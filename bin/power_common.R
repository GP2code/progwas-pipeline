#!/usr/bin/env Rscript

# Shared helpers for the power-analysis scripts (power_analysis_glm.R and friends).
# Base R only -- no optparse, no yaml, no ggplot2 -- so it is sourceable anywhere R runs.
# Sourced by the design scripts; not meant to be run directly.

# ==============================================================================
# Non-centrality-parameter core
#
# A 1-df association test is central chi-square under the null and non-central
# with parameter lambda = N_info * Var(G) * effect^2 under the alternative.
# These five functions are design-agnostic: only N_info differs between
# GLM (N), LMM (N*k/(1+(k-1)*rho)) and CoxPH (D events).
# ==============================================================================

# Lambda needed to reach target power. 39.601 at power=0.80, alpha=5e-8.
lambda_req <- function(power, alpha) (qnorm(1 - alpha / 2) + qnorm(power))^2

power_from_lambda <- function(lambda, alpha) {
  pchisq(qchisq(1 - alpha, df = 1), df = 1, ncp = lambda, lower.tail = FALSE)
}

# Genotypic variance, additive coding under HWE, scaled by imputation quality.
varG <- function(maf, info = 1) 2 * maf * (1 - maf) * info

# Minimum detectable variance explained. MAF-free, so comparable across designs.
mde_R2 <- function(N_info, alpha, power) lambda_req(power, alpha) / N_info

# MAF-dependent translation of mde_R2 onto the per-allele effect scale.
mde_beta <- function(N_info, maf, info, alpha, power) {
  sqrt(mde_R2(N_info, alpha, power) / varG(maf, info))
}

# ==============================================================================
# Study thresholds
# ==============================================================================

# Lowest MAF actually tested. At small N the MAC rule binds harder than the MAF rule.
maf_floor <- function(N_subjects, maf_thresh, mac_thresh) {
  max(maf_thresh, mac_thresh / (2 * N_subjects))
}

maf_floor_basis <- function(N_subjects, maf_thresh, mac_thresh) {
  sprintf("max(minor_allele_freq=%s, minor_allele_ct=%s/(2*%d)=%.4f)",
          format(maf_thresh), format(mac_thresh), N_subjects,
          mac_thresh / (2 * N_subjects))
}

# ==============================================================================
# Argument parsing helpers
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# Accepts "true"/"false"/"yes"/"no"/"1"/"0" so Nextflow can pass a boolean positively.
as_flag <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x)) return(default)
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "yes", "y", "1")
}

# Splits on whitespace and/or commas. Brackets and quotes are stripped first so a
# Groovy list interpolated by Nextflow ("[0.01, 0.05]") parses like "0.01 0.05".
parse_list <- function(x) {
  if (is.null(x) || !nzchar(trimws(as.character(x)))) return(character(0))
  cleaned <- gsub("[]['\"]", "", trimws(as.character(x)))
  parts <- unlist(strsplit(cleaned, "[[:space:],]+"))
  parts[nzchar(parts)]
}

parse_num_list <- function(x, label = "value list") {
  parts <- parse_list(x)
  if (!length(parts)) return(numeric(0))
  out <- suppressWarnings(as.numeric(parts))
  if (any(is.na(out))) {
    stop(sprintf("Could not parse %s as numbers: '%s'", label, paste(parts, collapse = " ")))
  }
  out
}

# ==============================================================================
# Input readers
# ==============================================================================

# Tab if the header has one, else comma. Matches the detection in main.nf.
detect_delim <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  header <- readLines(con, n = 1)
  if (!length(header) || !nzchar(header)) {
    stop(sprintf("'%s' appears to be empty (no header line found).", path))
  }
  if (grepl("\t", header, fixed = TRUE)) "\t" else ","
}

# Normalises PLINK-style #IID/#FID headers (read.table turns '#' into 'X.') and
# drops FID, since samples are keyed on IID alone. Mirrors bin/survival.R.
fix_plink_headers <- function(df, label) {
  names(df) <- sub("^X\\.IID$", "IID", names(df))
  names(df) <- sub("^X\\.FID$", "FID", names(df))
  if (!("IID" %in% names(df))) {
    stop(sprintf(paste0("'IID' column not found in --%s. Found columns: %s. ",
                        "The file requires an 'IID' column (or '#IID', which is also accepted)."),
                 label, paste(names(df), collapse = ", ")))
  }
  if ("FID" %in% names(df)) {
    df <- df[, names(df) != "FID", drop = FALSE]
  }
  df$IID <- trimws(as.character(df$IID))
  df
}

read_keyed_table <- function(path, label) {
  if (is.null(path) || !nzchar(path)) stop(sprintf("--%s is required.", label))
  if (!file.exists(path)) stop(sprintf("--%s not found: %s", label, path))
  # fill/na.strings tolerate ragged rows and empty trailing fields, which the
  # bundled example/phenotype.cs2.tsv has and which pandas accepts elsewhere in the pipeline.
  df <- read.table(path, header = TRUE, sep = detect_delim(path),
                   comment.char = "", check.names = TRUE, stringsAsFactors = FALSE,
                   fill = TRUE, na.strings = c("NA", "NaN", ""))
  fix_plink_headers(df, label)
}

# Merges on IID, failing if the covariate file has duplicate IIDs. Without this
# guard a Cartesian join silently inflates the frame (see survival.R commit bd70e45).
merge_pheno_covar <- function(pheno, covar) {
  dup <- sum(duplicated(covar$IID))
  if (dup > 0) {
    stop(sprintf(paste0("--covar-file contains %d duplicate IID(s); it must have one row per subject. ",
                        "First duplicate: '%s'."),
                 dup, covar$IID[duplicated(covar$IID)][1]))
  }
  shared <- setdiff(intersect(names(pheno), names(covar)), "IID")
  if (length(shared)) {
    covar <- covar[, !(names(covar) %in% shared), drop = FALSE]
  }
  merged <- merge(covar, pheno, by = "IID")
  if (nrow(merged) == 0) {
    stop("No IIDs are shared between --pheno-file and --covar-file.")
  }
  if (nrow(merged) > nrow(pheno)) {
    stop(sprintf("Merged data has %d rows but the phenotype file has %d; check for duplicate IIDs.",
                 nrow(merged), nrow(pheno)))
  }
  merged
}

# Splits into an 'ALL' pseudo-arm plus one entry per level of arm_col.
split_by_arm <- function(df, arm_col, per_arm = TRUE) {
  out <- list(ALL = df)
  if (!per_arm || is.null(arm_col) || !nzchar(arm_col) || !(arm_col %in% names(df))) {
    return(out)
  }
  levels_present <- sort(unique(as.character(df[[arm_col]])))
  levels_present <- levels_present[nzchar(levels_present) & !is.na(levels_present)]
  for (lv in levels_present) out[[lv]] <- df[as.character(df[[arm_col]]) == lv, , drop = FALSE]
  out
}

# ==============================================================================
# Phenotype scale
# ==============================================================================

# SD used to convert MDE beta from SD units onto the raw phenotype scale.
# Rank-INT phenotypes are already N(0,1), so BETA is in SD units by construction.
pheno_sd <- function(df, ycol, num_covs = character(0), cat_covs = character(0),
                     rank_int = FALSE, override = NA_real_, residualize = TRUE) {
  if (isTRUE(rank_int)) {
    return(list(sd = 1, source = "1 by construction (rank-INT applied before GWAS)"))
  }
  if (!is.na(override)) {
    return(list(sd = override, source = "user-supplied --sd-pheno (measured)"))
  }
  y <- suppressWarnings(as.numeric(df[[ycol]]))
  covs <- c(intersect(num_covs, names(df)), intersect(cat_covs, names(df)))
  if (!residualize || !length(covs)) {
    return(list(sd = sd(y, na.rm = TRUE), source = "marginal SD (not residualised)"))
  }
  fit_df <- df
  fit_df[[ycol]] <- y
  for (cc in intersect(cat_covs, names(fit_df))) fit_df[[cc]] <- factor(fit_df[[cc]])
  fml <- as.formula(paste(ycol, "~", paste(covs, collapse = " + ")))
  fit <- tryCatch(lm(fml, data = fit_df), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(sd = sd(y, na.rm = TRUE),
                source = "marginal SD (covariate model failed to fit)"))
  }
  list(sd = summary(fit)$sigma,
       source = paste0("residual SD after ", paste(covs, collapse = " + ")))
}

# ==============================================================================
# Output writers
#
# Three artefacts per run: the headline results table, an extended-info table
# holding everything that does not directly inform the headline, and a plot.
# One writer for every design, so GLM/LMM/CoxPH stay schema-compatible.
# ==============================================================================

# Expected minor allele count at a given MAF: 2N alleles per diploid sample.
mac_at_maf <- function(N_subjects, maf) 2 * N_subjects * maf

# One entry per group (study arm, or a single planned N). Designs fill what
# applies and leave the rest NA.
power_row <- function(Phenotype, Design, Group, N_info,
                      N_info_formula, N_info_basis,
                      Alpha, Target_power,
                      N_subjects = NA_integer_, N_obs = NA_integer_,
                      N_events = NA_integer_, Censoring_pct = NA_real_,
                      Visits_mean_k = NA_real_, ICC = NA_real_, ICC_source = NA_character_,
                      Info_r2 = 1, Info_source = NA_character_,
                      MAF_threshold = NA_real_, MAC_threshold = NA_real_,
                      MAF_floor = NA_real_, MAF_floor_basis = NA_character_,
                      SD_pheno = NA_real_, SD_source = NA_character_,
                      Rank_INT = NA, Effect_scale = NA_character_,
                      Model_package = NA_character_, Power_packages = NA_character_,
                      Power_method = NA_character_, Method_citation = NA_character_) {
  row <- as.list(environment())
  row$Lambda_req <- lambda_req(Target_power, Alpha)
  row$MDE_R2 <- row$Lambda_req / N_info
  row$MDE_R2_pct <- 100 * row$MDE_R2
  row
}

round_numeric <- function(df, digits = 6) {
  for (nm in names(df)) if (is.numeric(df[[nm]])) df[[nm]] <- signif(df[[nm]], digits)
  df
}

write_tsv <- function(df, path) {
  write.table(round_numeric(df), path, sep = "\t", row.names = FALSE,
              quote = FALSE, na = "NA")
  path
}

# Column glossaries. Written from the same name lists the tables are built from,
# so a new column cannot silently appear without a description.
COLUMN_KEY <- c(
  Phenotype        = "Phenotype column analysed; 'planned' in planning mode.",
  Design           = "Analysis design: glm, lmm_gallop or cph.",
  Group            = "Study arm, 'ALL' for the pooled cohort, or 'N=<n>' in planning mode.",
  N_info           = "Effective information: samples for GLM/LMM, observed events for CoxPH. Drives all MDEs.",
  MAF              = "Minor allele frequency evaluated at this row.",
  MAC              = "Expected minor allele count at this MAF: 2 * N * MAF.",
  MDE_R2_pct       = "Smallest variance explained (%) detectable at the target power. MAF-free, so constant within a group.",
  MDE_beta         = "Smallest per-allele effect detectable, in SD units of the phenotype.",
  MDE_beta_raw     = "MDE_beta rescaled to the phenotype's own units (MDE_beta * SD_pheno). NA if no SD is available.",
  MDE_HR           = "Smallest per-allele hazard ratio detectable (CoxPH only; NA otherwise).",
  Passes_QC        = "TRUE when this MAF survives the study filters: MAF >= MAF_threshold and MAC >= MAC_threshold.",
  N_subjects       = "Distinct subjects contributing to this group.",
  N_obs            = "Rows contributing (equals N_subjects for cross-sectional data).",
  N_events         = "Observed events (CoxPH only). This, not N_subjects, is what powers a survival analysis.",
  Censoring_pct    = "Percentage of subjects censored (CoxPH only).",
  Visits_mean_k    = "Mean repeated measurements per subject (LMM only).",
  ICC              = "Intraclass correlation of repeated measures (LMM only); higher ICC means repeat visits add less information.",
  ICC_source       = "How the ICC was obtained: estimated from the data or user-supplied.",
  N_info_formula   = "Formula used for N_info: 'N', 'N*k/(1+(k-1)*rho)' or 'D (events)'.",
  N_info_basis     = "Plain-language statement of exactly what was counted, including QC caveats.",
  Alpha            = "Significance threshold (5e-8 = genome-wide).",
  Target_power     = "Power the MDEs are solved for (0.80 = 80%).",
  Lambda_req       = "Non-centrality parameter required at this alpha and power; 39.601 at 5e-8 / 0.80.",
  MDE_R2           = "MDE_R2_pct as a proportion rather than a percentage.",
  Info_r2          = "Imputation quality r2 applied to Var(G). 1.0 means genotyped or assumed perfect.",
  Info_source      = "Whether Info_r2 was measured, user-supplied, or assumed.",
  MAF_threshold    = "Study MAF filter (params.minor_allele_freq), used for Passes_QC.",
  MAC_threshold    = "Study MAC filter (params.minor_allele_ct), used for Passes_QC.",
  MAF_floor        = "Lowest MAF actually testable: max(MAF_threshold, MAC_threshold/(2N)).",
  MAF_floor_basis  = "The MAF_floor arithmetic written out.",
  SD_pheno         = "Phenotype SD used to convert MDE_beta to raw units. 1 if rank-INT transformed.",
  SD_source        = "How SD_pheno was obtained: residual after covariates, marginal, supplied, or rank-INT.",
  Rank_INT         = "TRUE if the phenotype was rank-inverse-normal transformed before the GWAS.",
  Effect_scale     = "Units of the MDE effect columns.",
  Model_package    = "Software that fits the association model in the pipeline.",
  Power_packages   = "Software that computed these power numbers (distinct from Model_package).",
  Power_method     = "Statistical framework and its assumptions.",
  Method_citation  = "References for the method and the additive-coding convention."
)

write_column_key <- function(cols, path) {
  desc <- unname(COLUMN_KEY[cols])
  desc[is.na(desc)] <- "(no description available)"
  write_tsv(data.frame(Column = cols, Description = desc, stringsAsFactors = FALSE), path)
}

# Headline table: one row per group x MAF. MAC and Passes_QC show which MAF
# values the study's own thresholds actually keep -- the knobs a researcher tunes.
write_results_table <- function(rows, maf_grid, out_prefix, hazard = FALSE) {
  out <- do.call(rbind, lapply(rows, function(r) {
    mac <- mac_at_maf(if (is.na(r$N_subjects)) r$N_info else r$N_subjects, maf_grid)
    beta <- mde_beta(r$N_info, maf_grid, r$Info_r2, r$Alpha, r$Target_power)
    data.frame(
      Phenotype = r$Phenotype, Design = r$Design, Group = r$Group,
      N_info = r$N_info, MAF = maf_grid, MAC = mac,
      MDE_R2_pct = r$MDE_R2_pct,
      MDE_beta = beta,
      MDE_beta_raw = if (is.na(r$SD_pheno)) NA_real_ else beta * r$SD_pheno,
      MDE_HR = if (hazard) exp(beta) else NA_real_,
      Passes_QC = (is.na(r$MAF_threshold) | maf_grid >= r$MAF_threshold) &
                  (is.na(r$MAC_threshold) | mac >= r$MAC_threshold),
      stringsAsFactors = FALSE)
  }))
  write_column_key(names(out), file.path(dirname(out_prefix), "power_results_key.tsv"))
  write_tsv(out, paste0(out_prefix, "_power_results.tsv"))
}

EXTENDED_COLS <- c(
  "Phenotype", "Design", "Group",
  "N_subjects", "N_obs", "N_events", "Censoring_pct", "Visits_mean_k",
  "ICC", "ICC_source",
  "N_info", "N_info_formula", "N_info_basis",
  "Alpha", "Target_power", "Lambda_req", "MDE_R2",
  "Info_r2", "Info_source",
  "MAF_threshold", "MAC_threshold", "MAF_floor", "MAF_floor_basis",
  "SD_pheno", "SD_source", "Rank_INT", "Effect_scale",
  "Model_package", "Power_packages", "Power_method", "Method_citation"
)

# Everything needed to interpret or reproduce the headline numbers, kept out of
# the results table so that stays readable.
write_extended_info <- function(rows, out_prefix) {
  out <- do.call(rbind, lapply(rows, function(r) {
    as.data.frame(r[EXTENDED_COLS], stringsAsFactors = FALSE)
  }))
  write_column_key(names(out), file.path(dirname(out_prefix), "extended_info_key.tsv"))
  write_tsv(out, paste0(out_prefix, "_power_extended_info.tsv"))
}

# MDE R2 = lambda_req / N_info, so this is a straight line on log-log axes.
# Each group is marked on it, which answers "how much more N would I need?".
plot_power_curve <- function(rows, out_prefix, alpha, target_power) {
  path <- paste0(out_prefix, "_power_curve.png")

  # Log-scale ticks at 1/2/5 x powers of ten, written out in full rather than
  # in scientific notation (sample sizes) and without trailing zeros (percentages).
  log_ticks <- function(lo, hi) {
    e <- floor(log10(lo)):ceiling(log10(hi))
    t <- as.vector(outer(c(1, 2, 5), 10^e))
    sort(unique(t[t >= lo & t <= hi]))
  }
  fmt_tick <- function(x) {
    format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE, big.mark = ",")
  }

  ok <- tryCatch({
    n_vals <- vapply(rows, function(r) r$N_info, numeric(1))
    labels <- vapply(rows, function(r) as.character(r$Group), character(1))
    design <- as.character(rows[[1]]$Design)
    lam <- lambda_req(target_power, alpha)

    # CoxPH power is driven by events, not people; every other design by samples.
    x_label <- if (identical(design, "cph")) "N_events" else "N_samples"

    lo <- 10^floor(log10(min(n_vals) / 3))
    hi <- 10^ceiling(log10(max(n_vals) * 3))
    curve_n <- exp(seq(log(lo), log(hi), length.out = 200))
    curve_y <- 100 * lam / curve_n

    png(path, width = 1000, height = 700, res = 120,
        type = if (capabilities("cairo")) "cairo" else "Xlib")
    on.exit(dev.off(), add = TRUE)
    par(mar = c(4.5, 4.8, 3.5, 1))
    plot(curve_n, curve_y, log = "xy", type = "l", lwd = 2, col = "#2c7bb6",
         axes = FALSE, xlab = x_label,
         ylab = expression("Minimum detectable " * R^2 * " (%)"),
         main = sprintf("Detectable effect vs sample size (alpha=%.3g, power=%.2f)",
                        alpha, target_power))
    xt <- log_ticks(min(curve_n), max(curve_n))
    yt <- log_ticks(min(curve_y), max(curve_y))
    axis(1, at = xt, labels = fmt_tick(xt))
    axis(2, at = yt, labels = fmt_tick(yt), las = 1)
    box()
    abline(v = xt, h = yt, col = "grey88", lty = 3)
    lines(curve_n, curve_y, lwd = 2, col = "#2c7bb6")
    points(n_vals, 100 * lam / n_vals, pch = 19, col = "#d7191c", cex = 1.2)
    # Alternate label side so groups with similar N do not overprint each other.
    side <- rep(c(4, 2), length.out = length(n_vals))[rank(n_vals, ties.method = "first")]
    text(n_vals, 100 * lam / n_vals, labels = labels, pos = side, cex = 0.75, col = "#d7191c")
    TRUE
  }, error = function(e) {
    message("Plot skipped: ", conditionMessage(e)); FALSE
  })
  if (ok) invisible(path) else invisible(NULL)
}

# Records the environment and exact invocation so a reviewer can audit the numbers.
write_provenance <- function(out_prefix, opt, inputs = list(), extra = character(0)) {
  path <- paste0(out_prefix, "_provenance.txt")
  pkg_version <- function(p) {
    tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
  }
  invocation <- paste(commandArgs(trailingOnly = FALSE), collapse = " ")
  if (nchar(invocation) > 2000) invocation <- paste0(substr(invocation, 1, 2000), " ...[truncated]")
  lines <- c(
    "# proGWAS power analysis provenance",
    paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("R: ", R.version.string),
    paste0("optparse: ", pkg_version("optparse")),
    paste0("platform: ", paste(Sys.info()[c("sysname", "release", "nodename")], collapse = " / ")),
    paste0("invocation: ", invocation),
    "",
    "## resolved parameters"
  )
  for (nm in sort(names(opt))) {
    v <- opt[[nm]]
    if (is.null(v)) v <- "NULL"
    lines <- c(lines, sprintf("%s: %s", nm, paste(as.character(v), collapse = " ")))
  }
  if (length(inputs)) {
    lines <- c(lines, "", "## inputs")
    for (nm in names(inputs)) {
      p <- inputs[[nm]]
      if (!is.null(p) && nzchar(p) && file.exists(p)) {
        lines <- c(lines, sprintf("%s: %s (%d bytes, mtime %s)", nm, p,
                                  as.numeric(file.size(p)),
                                  format(file.mtime(p), "%Y-%m-%d %H:%M:%S")))
      } else {
        lines <- c(lines, sprintf("%s: %s (not supplied)", nm, if (is.null(p)) "NULL" else p))
      }
    }
  }
  if (length(extra)) lines <- c(lines, "", "## notes", extra)
  writeLines(lines, path)
  invisible(path)
}
