#!/usr/bin/env Rscript

# Power analysis for time-to-event GWAS (Cox proportional hazards).
# Survival power is driven by the number of observed EVENTS, not the number of
# people enrolled: N_info = D. Censored subjects contribute follow-up but no event.
#
# Two modes:
#   data mode     counts events in the phenotype file
#   planning mode --n-grid as enrolled samples with --event-rate, or as events directly

suppressPackageStartupMessages(library("optparse"))

# Nextflow puts bin/ on PATH but runs each task in its own work dir, so source()
# needs an absolute path. --file= is set by Rscript; Sys.which is the PATH fallback.
power_lib_dir <- function() {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) {
    d <- dirname(normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE))
    if (nzchar(d) && dir.exists(d)) return(d)
  }
  w <- Sys.which("power_analysis_cph.R")
  if (nzchar(w)) return(dirname(normalizePath(w)))
  getwd()
}

option_list <- list(
  make_option(c("--n-grid"), type = "character", default = NULL,
              help = paste("Planning mode: space-separated sample sizes to evaluate.",
                           "No phenotype or covariate file is needed.")),
  make_option(c("--event-col"), type = "character", default = NULL,
              help = "0/1 event indicator column; defaults to --pheno-name"),
  make_option(c("--event-rate"), type = "double", default = NA_real_,
              help = paste("Planning mode: fraction of enrolled subjects expected to have an event,",
                           "so D = N * rate. Without it, --n-grid values are read as event counts")),
  make_option(c("--pheno-file"), type = "character", default = NULL,
              help = "Phenotype file (TSV/CSV; IID or #IID column required)"),
  make_option(c("--covar-file"), type = "character", default = NULL,
              help = "Covariate file; one row per subject"),
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

if (planning_mode) {
  if (any(n_grid < 2)) stop("--n-grid values must be at least 2.")
  n_grid <- sort(unique(round(n_grid)))
  if (!is.na(opt[["event-rate"]]) && (opt[["event-rate"]] <= 0 || opt[["event-rate"]] > 1)) {
    stop("--event-rate must be in (0, 1].")
  }
} else if (is.null(opt[["pheno-file"]]) || is.null(opt[["covar-file"]])) {
  stop("Supply either --n-grid (planning mode) or both --pheno-file and --covar-file.")
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

# PCs are produced later by COMPUTE_PCA and are absent from a raw covariate file.
num_covs <- num_covs[!grepl("^(?i)PC[0-9]+$", num_covs, perl = TRUE)]

# ---- read -------------------------------------------------------------------

if (!planning_mode) {

pheno <- read_keyed_table(opt[["pheno-file"]], "pheno-file")
covar <- read_keyed_table(opt[["covar-file"]], "covar-file")

ecol <- if (!is.null(opt[["event-col"]])) opt[["event-col"]] else ycol
if (!(ecol %in% names(pheno))) {
  stop(sprintf("Event column '%s' not found in --pheno-file. Found columns: %s",
               ecol, paste(names(pheno), collapse = ", ")))
}
ycol <- ecol

merged <- merge_pheno_covar(pheno, covar)
merged[[ycol]] <- suppressWarnings(as.numeric(merged[[ycol]]))
merged <- merged[!is.na(merged[[ycol]]), , drop = FALSE]
if (nrow(merged) == 0) {
  stop(sprintf("--pheno-name '%s' has no non-missing numeric values after merging.", ycol))
}

# A Cox event indicator must be 0/1; anything else means the wrong column was named.
distinct_vals <- sort(unique(merged[[ycol]]))
if (!all(distinct_vals %in% c(0, 1))) {
  stop(sprintf(paste0("Event column '%s' must be a 0/1 indicator, but contains: %s. ",
                      "Name the event indicator with --event-col, or use ",
                      "power_analysis_glm.R for a quantitative phenotype."),
               ycol, paste(utils::head(distinct_vals, 6), collapse = ", ")))
}

# Counting-process format gives several intervals per subject but at most one event,
# so collapse to one row per subject carrying whether that subject ever had an event.
n_dup <- sum(duplicated(merged$IID))
ev_per_iid <- tapply(merged[[ycol]], merged$IID, function(v) as.integer(any(v == 1)))
multi <- sum(tapply(merged[[ycol]], merged$IID, function(v) sum(v == 1)) > 1)
if (multi > 0) {
  warning(sprintf("%d subject(s) carry more than one event row; counted once each.", multi))
}
first <- merged[!duplicated(merged$IID), , drop = FALSE]
first[[ycol]] <- as.integer(ev_per_iid[as.character(first$IID)])
merged <- first

}  # end data-mode input handling

# ---- compute ----------------------------------------------------------------

r_prov <- sprintf("R %s.%s base stats (qnorm/qchisq/pchisq); optparse %s. No genpwr.",
                  R.version$major, R.version$minor,
                  tryCatch(as.character(packageVersion("optparse")),
                           error = function(e) "unknown"))

common <- list(
  Phenotype = if (planning_mode) "planned" else ycol,
  Design = "cph",
  N_info_formula = "D (events)",
  Alpha = opt$alpha, Target_power = opt[["target-power"]],
  Info_r2 = opt$info,
  Info_source = if (opt$info == 1) "assumed 1.0 (no genotype files read)" else "user-supplied --info",
  MAF_threshold = opt[["minor-allele-freq"]], MAC_threshold = opt[["minor-allele-ct"]],
  Effect_scale = "log hazard ratio per allele",
  Model_package = "bin/survival.R - R survival::coxph",
  Power_packages = r_prov,
  Power_method = "Non-centrality-parameter framework (Schoenfeld); additive coding; HWE; 1 df; events-based",
  Method_citation = "Schoenfeld 1983 (events-based survival power); Sham & Purcell 2014 (NCP); Nalls et al. 2019 (additive coding)"
)

rows <- list()

if (planning_mode) {
  rate <- opt[["event-rate"]]
  for (n in n_grid) {
    # Without a rate the grid is already event counts; with one it is enrolment.
    d_events <- if (is.na(rate)) n else max(1, round(n * rate))
    n_enrolled <- if (is.na(rate)) NA_integer_ else n
    grp <- if (is.na(rate)) sprintf("D=%d", d_events) else sprintf("N=%d (D=%d)", n, d_events)
    # The MAF/MAC floor is a genotyping property, so it keys on people when known.
    n_for_floor <- if (is.na(rate)) d_events else n
    rows[[as.character(n)]] <- do.call(power_row, c(common, list(
      Group = grp, N_info = d_events,
      N_subjects = n_for_floor, N_events = d_events,
      Censoring_pct = if (is.na(rate)) NA_real_ else 100 * (1 - d_events / n),
      N_info_basis = if (is.na(rate))
          "planned event count from --n-grid (hypothetical; no enrolment assumed)"
        else sprintf("planned: %d enrolled x %g event rate = %d events", n, rate, d_events),
      MAF_floor = maf_floor(n_for_floor, opt[["minor-allele-freq"]], opt[["minor-allele-ct"]]),
      MAF_floor_basis = maf_floor_basis(n_for_floor, opt[["minor-allele-freq"]], opt[["minor-allele-ct"]]),
      SD_pheno = NA_real_, SD_source = "not applicable (effects are on the log-hazard scale)",
      Rank_INT = FALSE)))
  }
} else {
  arms <- split_by_arm(merged, opt[["study-arm-col"]], per_arm)
  for (arm in names(arms)) {
    d <- arms[[arm]]
    n <- length(unique(d$IID))
    d_events <- sum(d[[ycol]] == 1, na.rm = TRUE)
    if (d_events < 2) {
      warning(sprintf("Study arm '%s' has %d event(s); skipped. Cox power depends on events, not enrolment.",
                      arm, d_events))
      next
    }
    rows[[arm]] <- do.call(power_row, c(common, list(
      Group = arm, N_info = d_events,
      N_subjects = n, N_obs = nrow(d), N_events = d_events,
      Censoring_pct = 100 * (1 - d_events / n),
      N_info_basis = sprintf("%d observed events among %d subjects (%.1f%% censored). Power depends on events, NOT enrolment (PRE-QC)",
                             d_events, n, 100 * (1 - d_events / n)),
      MAF_floor = maf_floor(n, opt[["minor-allele-freq"]], opt[["minor-allele-ct"]]),
      MAF_floor_basis = maf_floor_basis(n, opt[["minor-allele-freq"]], opt[["minor-allele-ct"]]),
      SD_pheno = NA_real_, SD_source = "not applicable (effects are on the log-hazard scale)",
      Rank_INT = FALSE)))
  }
}

if (!length(rows)) stop("No group had enough subjects to compute power.")

results_path  <- write_results_table(rows, maf_grid, opt[["out-prefix"]], hazard = TRUE)
extended_path <- write_extended_info(rows, opt[["out-prefix"]])
plot_path     <- if (as_flag(opt$plots, TRUE)) {
  plot_power_curve(rows, opt[["out-prefix"]], opt$alpha, opt[["target-power"]])
} else NULL

write_provenance(
  opt[["out-prefix"]], opt,
  inputs = list(`pheno-file` = opt[["pheno-file"]], `covar-file` = opt[["covar-file"]]),
  extra = c(
    sprintf("mode: %s", if (planning_mode) "planning (--n-grid)" else "data (phenotype + covariate files)"),
    sprintf("lambda_req(%.4g, %.3g) = %.4f", opt[["target-power"]], opt$alpha,
            lambda_req(opt[["target-power"]], opt$alpha)),
    sprintf("MAF grid: %s", paste(maf_grid, collapse = " ")),
    "MDE_R2 is MAF-free (lambda_req / N_info); only MDE_beta varies across MAF.",
    "MAC = 2 * N * MAF; Passes_QC is TRUE where MAF >= minor_allele_freq and MAC >= minor_allele_ct.",
    "N_info = D (observed events). Enrolled N affects only the MAF/MAC floor, not power.",
    if (planning_mode) "Planning mode: N values are hypothetical, not read from data."
    else c(sprintf("repeated IID rows dropped: %d", n_dup),
           "N is counted from the input files (PRE-QC); the analysed N after ancestry/kinship/",
           "missingness filtering will be lower, so these MDEs are optimistic.")
  ))

# ---- summary ----------------------------------------------------------------

cat(sprintf("\nCoxPH power analysis - %s mode (alpha=%.3g, target power=%.2f)\n",
            if (planning_mode) "planning" else "data", opt$alpha, opt[["target-power"]]))
cat(sprintf("%-20s %9s %12s %14s\n", "group", "D (events)", "MDE R2 (%)", "MDE HR@0.05"))
for (nm in names(rows)) {
  r <- rows[[nm]]
  cat(sprintf("%-20s %9d %12.4f %14.4f\n", r$Group, as.integer(r$N_info), r$MDE_R2_pct,
              exp(mde_beta(r$N_info, 0.05, r$Info_r2, r$Alpha, r$Target_power))))
}
cat(sprintf("\nWrote: %s\n       %s\n", results_path, extended_path))
if (!is.null(plot_path)) cat(sprintf("       %s\n", plot_path))
cat(sprintf("       %s_provenance.txt\n", opt[["out-prefix"]]))
