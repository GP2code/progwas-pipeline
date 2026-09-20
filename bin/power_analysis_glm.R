#!/usr/bin/env Rscript

# Power analysis for cross-sectional quantitative GWAS (plink2 --glm), where the
# information carried by the design is just the number of subjects: N_info = N.
# Reads the phenotype and covariate files only; genotypes are not used.

suppressPackageStartupMessages(library("optparse"))

# Nextflow puts bin/ on PATH but runs each task in its own work dir, so source()
# needs an absolute path. --file= is set by Rscript; Sys.which is the PATH fallback.
power_lib_dir <- function() {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) {
    d <- dirname(normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE))
    if (nzchar(d) && dir.exists(d)) return(d)
  }
  w <- Sys.which("power_analysis_glm.R")
  if (nzchar(w)) return(dirname(normalizePath(w)))
  getwd()
}

option_list <- list(
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

if (is.null(opt[["pheno-file"]]) || is.null(opt[["covar-file"]])) {
  stop("--pheno-file and --covar-file are both required.")
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

pheno <- read_keyed_table(opt[["pheno-file"]], "pheno-file")
covar <- read_keyed_table(opt[["covar-file"]], "covar-file")

if (!(ycol %in% names(pheno))) {
  stop(sprintf("--pheno-name '%s' not found in --pheno-file. Found columns: %s",
               ycol, paste(names(pheno), collapse = ", ")))
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

# One row per subject is expected; a longitudinal file would otherwise inflate N.
n_dup <- sum(duplicated(merged$IID))
if (n_dup > 0) {
  warning(sprintf(paste0("--pheno-file has %d repeated IID row(s); keeping the first row per ",
                         "subject. For repeated measures use power_analysis_lmm.R."), n_dup))
  merged <- merged[!duplicated(merged$IID), , drop = FALSE]
}

# ---- compute ----------------------------------------------------------------

r_prov <- sprintf("R %s.%s base stats (qnorm/qchisq/pchisq); optparse %s. No genpwr.",
                  R.version$major, R.version$minor,
                  tryCatch(as.character(packageVersion("optparse")),
                           error = function(e) "unknown"))

arms <- split_by_arm(merged, opt[["study-arm-col"]], per_arm)
rows <- list()

for (arm in names(arms)) {
  d <- arms[[arm]]
  n <- length(unique(d$IID))
  if (n < 2) {
    warning(sprintf("Study arm '%s' has %d subject(s); skipped.", arm, n))
    next
  }
  sd_info <- pheno_sd(d, ycol, num_covs, cat_covs,
                      rank_int = opt[["rank-int"]], override = opt[["sd-pheno"]],
                      residualize = residualize)
  rows[[arm]] <- power_row(
    Phenotype = ycol, Design = "glm", Term = "main", Study_arm = arm,
    N_subjects = n, N_obs = nrow(d),
    N_info = n, N_info_formula = "N",
    N_info_basis = "unique IIDs with non-missing phenotype, present in the covariate file (PRE-QC)",
    Alpha = opt$alpha, Target_power = opt[["target-power"]],
    Info_r2 = opt$info,
    Info_source = if (opt$info == 1) "assumed 1.0 (no genotype files read)" else "user-supplied --info",
    MAF_floor = maf_floor(n, opt[["minor-allele-freq"]], opt[["minor-allele-ct"]]),
    MAF_floor_basis = maf_floor_basis(n, opt[["minor-allele-freq"]], opt[["minor-allele-ct"]]),
    SD_pheno = sd_info$sd, SD_source = sd_info$source,
    Rank_INT = opt[["rank-int"]], Effect_scale = "SD units of phenotype per allele",
    Model_package = "plink2 --glm (additive, 1 df)",
    Power_packages = r_prov,
    Power_method = "Non-centrality-parameter framework; additive coding; HWE; 1 df Wald/chi-square",
    Method_citation = "Sham & Purcell 2014 (NCP); Nalls et al. 2019 (additive coding convention)"
  )
}

if (!length(rows)) stop("No study arm had enough subjects to compute power.")

paths <- write_power_tables(rows, maf_grid, opt[["out-prefix"]], hazard = FALSE)

write_provenance(
  opt[["out-prefix"]], opt,
  inputs = list(`pheno-file` = opt[["pheno-file"]], `covar-file` = opt[["covar-file"]]),
  extra = c(
    sprintf("lambda_req(%.4g, %.3g) = %.4f", opt[["target-power"]], opt$alpha,
            lambda_req(opt[["target-power"]], opt$alpha)),
    sprintf("MAF grid: %s", paste(maf_grid, collapse = " ")),
    sprintf("repeated IID rows dropped: %d", n_dup),
    "N is counted from the input files (PRE-QC); the analysed N after ancestry/kinship/",
    "missingness filtering will be lower, so these MDEs are optimistic."
  ))

# ---- summary ----------------------------------------------------------------

cat(sprintf("\nGLM power analysis: %s (alpha=%.3g, target power=%.2f)\n",
            ycol, opt$alpha, opt[["target-power"]]))
cat(sprintf("%-16s %8s %12s %14s\n", "study_arm", "N", "MDE R2 (%)", "MDE beta@0.05"))
for (arm in names(rows)) {
  r <- rows[[arm]]
  cat(sprintf("%-16s %8d %12.4f %14.4f\n", arm, r$N_subjects, r$MDE_R2_pct,
              mde_beta(r$N_info, 0.05, r$Info_r2, r$Alpha, r$Target_power)))
}
cat(sprintf("\nWrote: %s\n       %s\n       %s_provenance.txt\n",
            paths$wide, paths$long, opt[["out-prefix"]]))
