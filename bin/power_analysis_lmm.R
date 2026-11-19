#!/usr/bin/env Rscript

# Power analysis for longitudinal quantitative GWAS (GALLOP / linear mixed model).
# Repeat visits on one subject are correlated, so k visits are worth less than k
# independent observations: N_info = N*k / (1 + (k-1)*rho), with rho the ICC.
#
# Two modes:
#   data mode     reads phenotype + covariate files; k and rho estimated from the data
#   planning mode --n-grid with --visits and --icc, for sizing a study in advance

suppressPackageStartupMessages(library("optparse"))

# Nextflow puts bin/ on PATH but runs each task in its own work dir, so source()
# needs an absolute path. --file= is set by Rscript; Sys.which is the PATH fallback.
power_lib_dir <- function() {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) {
    d <- dirname(normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE))
    if (nzchar(d) && dir.exists(d)) return(d)
  }
  w <- Sys.which("power_analysis_lmm.R")
  if (nzchar(w)) return(dirname(normalizePath(w)))
  getwd()
}

option_list <- list(
  make_option(c("--n-grid"), type = "character", default = NULL,
              help = paste("Planning mode: space-separated sample sizes to evaluate.",
                           "No phenotype or covariate file is needed.")),
  make_option(c("--time-col"), type = "character", default = "study_days",
              help = "Visit/time column in the phenotype file [default %default]"),
  make_option(c("--icc"), type = "double", default = NA_real_,
              help = "Intraclass correlation (rho). Estimated from the data if omitted; required in planning mode"),
  make_option(c("--visits"), type = "double", default = NA_real_,
              help = "Planning mode: mean visits per subject (k)"),
  make_option(c("--min-visits"), type = "integer", default = 1,
              help = paste("Drop subjects with fewer rows than this before estimating the ICC.",
                           "The default keeps everyone: the one-way ANOVA already handles",
                           "single-visit subjects correctly (they inform the between-subject",
                           "variance and carry zero within-subject df) [default %default]")),
  make_option(c("--slope-term"), type = "character", default = "true",
              help = "Also report the SNPxtime (slope) term as an upper bound [default %default]"),
  make_option(c("--pheno-file"), type = "character", default = NULL,
              help = "Phenotype file (TSV/CSV; IID or #IID column required)"),
  make_option(c("--covar-file"), type = "character", default = NULL,
              help = "Covariate file; one row per subject"),
  make_option(c("--analysis-set"), type = "character", default = NULL,
              help = paste("Post-QC analysis set from a completed run:",
                           "prepared_data/*_filtered.pca.harmonized.tsv, comma-separated for",
                           "several study arms. Used INSTEAD of --covar-file, so N and the",
                           "covariates (PCs included) are the ones the GWAS actually used")),
  make_option(c("--pheno-name"), type = "character", default = "y",
              help = "Phenotype column to analyse [default %default]"),
  make_option(c("--study-arm-col"), type = "character", default = "study_arm",
              help = "Covariate column defining arms analysed separately [default %default]"),
  make_option(c("--covar-numeric"), type = "character", default = "",
              help = "Space-delimited numeric covariates, used for the residual SD"),
  make_option(c("--covar-categorical"), type = "character", default = "",
              help = "Space-delimited categorical covariates, used for the residual SD"),
  make_option(c("--maf-grid"), type = "character", default = "0.01 0.05 0.10 0.20 0.50",
              help = "MAF grid for the MDE columns [default %default]"),
  make_option(c("--minor-allele-freq"), type = "double", default = 0.05,
              help = "Study MAF threshold, for the effective MAF floor [default %default]"),
  make_option(c("--minor-allele-ct"), type = "double", default = 20,
              help = "Study MAC threshold, for the effective MAF floor [default %default]"),
  make_option(c("--alpha"), type = "double", default = 5e-8,
              help = "Significance threshold [default %default]"),
  make_option(c("--target-power"), type = "double", default = 0.80,
              help = "Target power [default %default]"),
  make_option(c("--info"), type = "double", default = 1.0,
              help = "Imputation r2 applied to Var(G); 1.0 = genotyped/perfectly imputed [default %default]"),
  make_option(c("--sd-pheno"), type = "double", default = NA_real_,
              help = "Measured phenotype SD; overrides the value computed from the data"),
  make_option(c("--rank-int"), action = "store_true", default = FALSE,
              help = "Phenotype was rank-INT transformed before GWAS, so BETA is already in SD units"),
  make_option(c("--residualize"), type = "character", default = "true",
              help = "Regress covariates out before taking the SD [default %default]"),
  make_option(c("--plots"), type = "character", default = "true",
              help = "Write the power-curve PNG [default %default]"),
  make_option(c("--per-arm"), type = "character", default = "true",
              help = "Also report each study arm separately [default %default]"),
  make_option(c("--out-prefix"), type = "character", default = "POWER",
              help = "Output file prefix [default %default]"),
  make_option(c("--power-lib"), type = "character", default = NULL,
              help = "Explicit path to power_common.R (normally self-located)")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = paste0(
    "\nMinimum detectable effect sizes for a cross-sectional quantitative GWAS.\n",
    "Power is computed analytically from the non-centrality parameter\n",
    "lambda = N * Var(G) * beta^2, assuming additive coding and HWE.\n")))

lib <- opt[["power-lib"]]
if (is.null(lib)) lib <- file.path(power_lib_dir(), "power_common.R")
if (!file.exists(lib)) {
  stop(sprintf("power_common.R not found (looked at: %s). Pass --power-lib explicitly.", lib))
}
source(lib)

# ---- validate ---------------------------------------------------------------

n_grid <- parse_num_list(opt[["n-grid"]], "--n-grid")
planning_mode <- length(n_grid) > 0

# A post-QC analysis set stands in for the raw covariate file: it holds exactly
# the samples the GWAS analysed, with the PCs COMPUTE_PCA produced merged in.
analysis_set <- parse_paths(opt[["analysis-set"]])
use_analysis_set <- length(analysis_set) > 0 && !planning_mode
if (use_analysis_set && !is.null(opt[["covar-file"]])) {
  warning("Both --analysis-set and --covar-file given; using --analysis-set (post-QC).")
}

if (planning_mode) {
  if (any(n_grid < 2)) stop("--n-grid values must be at least 2.")
  n_grid <- sort(unique(round(n_grid)))
  if (is.na(opt$visits) || opt$visits < 1) {
    stop("Planning mode needs --visits (mean measurements per subject); k=1 is a cross-sectional design, use power_analysis_glm.R.")
  }
  if (is.na(opt$icc) || opt$icc < 0 || opt$icc >= 1) {
    stop("Planning mode needs --icc in [0, 1): the ICC cannot be estimated without data. Typical repeated clinical measures fall around 0.3-0.6.")
  }
} else if (is.null(opt[["pheno-file"]]) || (is.null(opt[["covar-file"]]) && !use_analysis_set)) {
  stop(paste("Supply either --n-grid (planning mode), or --pheno-file plus one of",
             "--covar-file (pre-QC) and --analysis-set (post-QC, recommended after a run)."))
}
if (opt$alpha <= 0 || opt$alpha >= 1) stop("--alpha must be between 0 and 1.")
if (opt[["target-power"]] <= 0 || opt[["target-power"]] >= 1) {
  stop("--target-power must be between 0 and 1.")
}
if (opt$info <= 0 || opt$info > 1) stop("--info must be in (0, 1].")

maf_grid <- parse_num_list(opt[["maf-grid"]], "--maf-grid")
if (!length(maf_grid)) stop("--maf-grid is empty.")
if (any(maf_grid <= 0 | maf_grid >= 0.5 + 1e-9)) {
  stop("--maf-grid values must be in (0, 0.5].")
}
maf_grid <- sort(unique(maf_grid))

ycol <- opt[["pheno-name"]]
residualize <- as_flag(opt$residualize, TRUE)
per_arm <- as_flag(opt[["per-arm"]], TRUE)
num_covs <- parse_list(opt[["covar-numeric"]])
cat_covs <- parse_list(opt[["covar-categorical"]])

# PCs are absent from a raw covariate file (COMPUTE_PCA makes them later), but a
# post-QC analysis set already carries them, so keep them there: the residual SD
# should be taken over the same covariates the GWAS adjusted for.
if (!use_analysis_set) {
  num_covs <- num_covs[!grepl("^(?i)PC[0-9]+$", num_covs, perl = TRUE)]
}

# ---- read -------------------------------------------------------------------

if (!planning_mode) {

pheno <- read_keyed_table(opt[["pheno-file"]], "pheno-file")
covar <- if (use_analysis_set) read_analysis_set(analysis_set)
         else read_keyed_table(opt[["covar-file"]], "covar-file")

if (!(ycol %in% names(pheno))) {
  stop(sprintf("--pheno-name '%s' not found in --pheno-file. Found columns: %s",
               ycol, paste(names(pheno), collapse = ", ")))
}

tcol <- opt[["time-col"]]
if (!(tcol %in% names(pheno))) {
  stop(sprintf("--time-col '%s' not found in --pheno-file. Found columns: %s",
               tcol, paste(names(pheno), collapse = ", ")))
}

merged <- merge_pheno_covar(pheno, covar)
merged[[ycol]] <- suppressWarnings(as.numeric(merged[[ycol]]))
merged <- merged[!is.na(merged[[ycol]]), , drop = FALSE]
if (nrow(merged) == 0) {
  stop(sprintf("--pheno-name '%s' has no non-missing numeric values after merging.", ycol))
}

# Binary phenotypes need the case/control information formula, not N.
distinct_vals <- unique(merged[[ycol]])
if (length(distinct_vals) <= 2 && all(distinct_vals %in% c(0, 1))) {
  stop(sprintf(paste0("--pheno-name '%s' looks binary (only 0/1 values). Binary/logistic power ",
                      "is not supported here: N_info for a case/control contrast is ",
                      "n1*n2/(n1+n2), not N. Use a quantitative phenotype."), ycol))
}

# Repeated rows are the point here; without them there is no within-subject information.
n_dup <- sum(duplicated(merged$IID))
if (n_dup == 0) {
  stop(paste0("--pheno-file has one row per subject, so there are no repeated measures to model. ",
              "For cross-sectional data use power_analysis_glm.R."))
}

}  # end data-mode input handling

# ---- compute ----------------------------------------------------------------

r_prov <- sprintf("R %s.%s base stats (qnorm/qchisq/pchisq); optparse %s. No genpwr.",
                  R.version$major, R.version$minor,
                  tryCatch(as.character(packageVersion("optparse")),
                           error = function(e) "unknown"))

qc_tag <- if (use_analysis_set) "POST-QC" else "PRE-QC"

mode_label <- if (planning_mode) {
  "planning (--n-grid)"
} else if (use_analysis_set) {
  "after a run (post-QC analysis set + phenotype file)"
} else {
  "data (phenotype + covariate files, pre-QC)"
}

sample_source <- if (planning_mode) {
  "planning mode (no samples read)"
} else if (use_analysis_set) {
  sprintf("POST-QC analysis set from a completed run: %s",
          paste(basename(analysis_set), collapse = ", "))
} else {
  "PRE-QC covariate file (ancestry, kinship and missingness filters not yet applied)"
}

common <- list(
  Phenotype = if (planning_mode) "planned" else ycol,
  Design = "lmm_gallop",
  N_info_formula = "N*k/(1+(k-1)*rho)",
  Sample_source = sample_source,
  Alpha = opt$alpha, Target_power = opt[["target-power"]],
  Info_r2 = opt$info,
  Info_source = if (opt$info == 1) "assumed 1.0 (no genotype files read)" else "user-supplied --info",
  MAF_threshold = opt[["minor-allele-freq"]], MAC_threshold = opt[["minor-allele-ct"]],
  Effect_scale = "SD units of phenotype per allele",
  Model_package = "bin/gallop.py - GALLOP (Python: numpy, scipy, scikit-learn)",
  Power_packages = r_prov,
  Power_method = "Non-centrality-parameter framework; additive coding; HWE; 1 df; repeated-measures design effect",
  Method_citation = "Sham & Purcell 2014 (NCP); Nalls et al. 2019 (additive coding convention)"
)

slope_term <- as_flag(opt[["slope-term"]], TRUE)

# GALLOP's slope term is a different estimand from the main term. With
# subject-specific random slopes it carries at most one independent observation
# per subject, so N is an upper bound -- honest, and better than reusing the
# main-term number. The true value needs a fitted mixed model.
add_slope_row <- function(rows, key, group, n_subj, base) {
  if (!slope_term) return(rows)
  rows[[paste0(key, ":slope")]] <- do.call(power_row, c(common, base, list(
    Group = paste0(group, ":slope"), N_info = n_subj, N_subjects = n_subj,
    N_info_basis = paste("N (UPPER BOUND for the SNPxtime slope term, the interaction:",
                         "BETA_INT/SE_INT/P_INT in GALLOP output. True slope information is lower and",
                         "depends on visit-time spread and random-slope variance, which base R cannot",
                         "estimate without fitting a mixed model)"))))
  rows
}

rows <- list()

if (planning_mode) {
  k <- opt$visits
  rho <- opt$icc
  for (n in n_grid) {
    ninfo <- n_info_lmm(n, k, rho)
    base <- list(
      MAF_floor = maf_floor(n, opt[["minor-allele-freq"]], opt[["minor-allele-ct"]]),
      MAF_floor_basis = maf_floor_basis(n, opt[["minor-allele-freq"]], opt[["minor-allele-ct"]]),
      SD_pheno = opt[["sd-pheno"]],
      SD_source = if (is.na(opt[["sd-pheno"]])) "not available in planning mode" else "user-supplied --sd-pheno",
      Rank_INT = opt[["rank-int"]], Visits_mean_k = k, ICC = rho,
      ICC_source = "user-supplied --icc (planning mode)")
    rows[[as.character(n)]] <- do.call(power_row, c(common, base, list(
      Group = sprintf("N=%d:intercept", n), N_info = ninfo, N_subjects = n, N_obs = n * k,
      N_info_basis = sprintf(paste("planned: %d subjects x %g visits, deflated by ICC %g.",
                                   "Main term = the intercept (BETA/SE/P in GALLOP output)"), n, k, rho))))
    rows <- add_slope_row(rows, as.character(n), sprintf("N=%d", n), n, base)
  }
} else {
  arms <- split_by_arm(merged, opt[["study-arm-col"]], per_arm)
  for (arm in names(arms)) {
    d <- arms[[arm]]
    n <- length(unique(d$IID))
    if (n < 2) {
      warning(sprintf("Study arm '%s' has %d subject(s); skipped.", arm, n))
      next
    }
    # Subjects with too few visits carry no within-subject information for the ICC.
    counts <- table(d$IID)
    keep <- names(counts)[counts >= opt[["min-visits"]]]
    d_icc <- d[d$IID %in% keep, , drop = FALSE]
    est <- icc_oneway(d_icc[[ycol]], d_icc$IID)

    rho <- if (!is.na(opt$icc)) opt$icc else est$icc
    if (is.na(rho)) {
      warning(sprintf("Study arm '%s': ICC not estimable; skipped. Supply --icc to override.", arm))
      next
    }
    k <- length(d[[ycol]]) / n
    ninfo <- n_info_lmm(n, k, rho)

    # For repeated measures the beta scale is the observation-level total SD,
    # sqrt(var_between + var_within), not sd(y) of the stacked rows.
    sd_val <- if (!is.na(opt[["sd-pheno"]])) opt[["sd-pheno"]]
              else if (opt[["rank-int"]]) 1
              else if (!is.na(est$sd_total)) est$sd_total
              else sd(d[[ycol]], na.rm = TRUE)
    sd_src <- if (!is.na(opt[["sd-pheno"]])) "user-supplied --sd-pheno (measured)"
              else if (opt[["rank-int"]]) "1 by construction (rank-INT applied before GWAS)"
              else "sqrt(var_between + var_within) from the variance components"

    base <- list(
      MAF_floor = maf_floor(n, opt[["minor-allele-freq"]], opt[["minor-allele-ct"]]),
      MAF_floor_basis = maf_floor_basis(n, opt[["minor-allele-freq"]], opt[["minor-allele-ct"]]),
      SD_pheno = sd_val, SD_source = sd_src, Rank_INT = opt[["rank-int"]],
      Visits_mean_k = k, ICC = rho,
      ICC_source = if (!is.na(opt$icc)) "user-supplied --icc"
                   else sprintf("one-way random-effects ANOVA on the raw phenotype (%d subjects with >=%d visits)%s",
                                est$N, opt[["min-visits"]],
                                if (nzchar(est$note)) paste0("; ", est$note) else ""))
    rows[[arm]] <- do.call(power_row, c(common, base, list(
      Group = paste0(arm, ":intercept"), N_info = ninfo, N_subjects = n, N_obs = nrow(d),
      N_info_basis = sprintf(paste("%d subjects x %.3f mean visits, deflated by ICC %.4f (%s).",
                                   "Main term = the intercept (BETA/SE/P in GALLOP output)"),
                             n, k, rho, qc_tag))))
    rows <- add_slope_row(rows, arm, arm, n, base)
  }
}

if (!length(rows)) stop("No group had enough subjects to compute power.")

results_path  <- write_results_table(rows, maf_grid, opt[["out-prefix"]], hazard = FALSE)
extended_path <- write_extended_info(rows, opt[["out-prefix"]])
plot_path     <- if (as_flag(opt$plots, TRUE)) {
  plot_power_curve(rows, opt[["out-prefix"]], opt$alpha, opt[["target-power"]])
} else NULL

write_provenance(
  opt[["out-prefix"]], opt,
  inputs = c(list(`pheno-file` = opt[["pheno-file"]]),
             if (use_analysis_set)
               setNames(as.list(analysis_set), sprintf("analysis-set[%d]", seq_along(analysis_set)))
             else list(`covar-file` = opt[["covar-file"]])),
  extra = c(
    sprintf("mode: %s", mode_label),
    sprintf("lambda_req(%.4g, %.3g) = %.4f", opt[["target-power"]], opt$alpha,
            lambda_req(opt[["target-power"]], opt$alpha)),
    sprintf("MAF grid: %s", paste(maf_grid, collapse = " ")),
    "N_info = N*k/(1+(k-1)*rho): repeat visits on one subject are correlated, so k visits",
    "are worth less than k independent observations.",
    if (as_flag(opt[["slope-term"]], TRUE))
      "':slope' rows bound the SNPxtime term at N; the true value needs a fitted mixed model."
    else NULL,
    "MDE_R2 is MAF-free (lambda_req / N_info); only MDE_beta varies across MAF.",
    "MAC = 2 * N * MAF; Passes_QC is TRUE where MAF >= minor_allele_freq and MAC >= minor_allele_ct.",
    if (planning_mode) "Planning mode: N values are hypothetical, not read from data."
    else c(sprintf("repeated IID rows (the point of a longitudinal design): %d", n_dup),
           if (use_analysis_set)
             c("N comes from the post-QC analysis set, so it is the number of samples the GWAS",
               "actually analysed, and the covariates include the PCs the pipeline computed.")
           else
             c("N is counted from the input files (PRE-QC); the analysed N after ancestry/kinship/",
               "missingness filtering will be lower, so these MDEs are optimistic."))
  ))

# ---- summary ----------------------------------------------------------------

cat(sprintf("\nLMM (GALLOP) power analysis - %s mode (alpha=%.3g, target power=%.2f)\n",
            if (planning_mode) "planning" else if (use_analysis_set) "after-run" else "data",
            opt$alpha, opt[["target-power"]]))
cat(sprintf("%-24s %9s %12s %14s\n", "group", "N_info", "MDE R2 (%)", "MDE beta@0.05"))
for (nm in names(rows)) {
  r <- rows[[nm]]
  cat(sprintf("%-24s %9d %12.4f %14.4f\n", r$Group, as.integer(r$N_info), r$MDE_R2_pct,
              mde_beta(r$N_info, 0.05, r$Info_r2, r$Alpha, r$Target_power)))
}
cat(sprintf("\nWrote: %s\n       %s\n", results_path, extended_path))
if (!is.null(plot_path)) cat(sprintf("       %s\n", plot_path))
cat(sprintf("       %s_provenance.txt\n", opt[["out-prefix"]]))
