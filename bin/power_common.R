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
  df <- read.table(path, header = TRUE, sep = detect_delim(path),
                   comment.char = "", check.names = TRUE, stringsAsFactors = FALSE)
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
# One writer for every design, so GLM/LMM/CoxPH rows land in a single table with
# identical columns, ordering and rounding.
# ==============================================================================

# Fixed scalar schema. Designs fill what applies and leave the rest NA.
power_row <- function(Phenotype, Design, Term = "main", Study_arm = "ALL",
                      N_subjects = NA_integer_, N_obs = NA_integer_,
                      N_events = NA_integer_, Censoring_pct = NA_real_,
                      Visits_mean_k = NA_real_, ICC = NA_real_, ICC_source = NA_character_,
                      N_info = NA_real_, N_info_formula = NA_character_,
                      N_info_basis = NA_character_,
                      Alpha = NA_real_, Target_power = NA_real_,
                      Info_r2 = NA_real_, Info_source = NA_character_,
                      MAF_floor = NA_real_, MAF_floor_basis = NA_character_,
                      SD_pheno = NA_real_, SD_source = NA_character_,
                      Rank_INT = NA, Effect_scale = NA_character_,
                      Model_package = NA_character_, Power_packages = NA_character_,
                      Power_method = NA_character_, Method_citation = NA_character_) {
  row <- as.list(environment())
  row$Lambda_req <- if (!is.na(Alpha) && !is.na(Target_power)) lambda_req(Target_power, Alpha) else NA_real_
  row$MDE_R2 <- if (!is.na(N_info) && !is.na(row$Lambda_req)) row$Lambda_req / N_info else NA_real_
  row$MDE_R2_pct <- 100 * row$MDE_R2
  row
}

POWER_SCALAR_COLS <- c(
  "Phenotype", "Design", "Term", "Study_arm",
  "N_subjects", "N_obs", "N_events", "Censoring_pct", "Visits_mean_k",
  "ICC", "ICC_source",
  "N_info", "N_info_formula", "N_info_basis",
  "Alpha", "Target_power", "Lambda_req", "MDE_R2", "MDE_R2_pct",
  "Info_r2", "Info_source", "MAF_floor", "MAF_floor_basis",
  "SD_pheno", "SD_source", "Rank_INT", "Effect_scale",
  "Model_package", "Power_packages", "Power_method", "Method_citation"
)

round_numeric <- function(df, digits = 6) {
  for (nm in names(df)) {
    if (is.numeric(df[[nm]])) df[[nm]] <- signif(df[[nm]], digits)
  }
  df
}

# Column-name suffix for a MAF. %g avoids the leading-space padding formatC adds,
# which would otherwise split the TSV header.
maf_label <- function(maf) gsub("[^0-9.]", "", sprintf("%g", maf))

# Writes the wide per-design table (one row per arm, MDE columns across the grid)
# and the tidy long table (one row per arm x MAF).
write_power_tables <- function(rows, maf_grid, out_prefix, hazard = FALSE) {
  stopifnot(length(rows) > 0)

  scalars <- do.call(rbind, lapply(rows, function(r) {
    as.data.frame(r[POWER_SCALAR_COLS], stringsAsFactors = FALSE)
  }))

  wide <- scalars
  long <- list()

  for (i in seq_along(rows)) {
    r <- rows[[i]]
    info <- if (is.na(r$Info_r2)) 1 else r$Info_r2
    beta <- mde_beta(r$N_info, maf_grid, info, r$Alpha, r$Target_power)
    beta_raw <- if (is.na(r$SD_pheno)) rep(NA_real_, length(maf_grid)) else beta * r$SD_pheno
    hr <- if (hazard) exp(beta) else rep(NA_real_, length(maf_grid))

    long[[i]] <- data.frame(
      Phenotype = r$Phenotype, Design = r$Design, Term = r$Term, Study_arm = r$Study_arm,
      N_info = r$N_info, MAF = maf_grid, Var_G = varG(maf_grid, info), Info_r2 = info,
      MDE_R2 = r$MDE_R2, MDE_beta = beta, MDE_beta_raw = beta_raw, MDE_HR = hr,
      Below_MAF_floor = maf_grid < r$MAF_floor,
      Alpha = r$Alpha, Target_power = r$Target_power,
      stringsAsFactors = FALSE
    )

    for (j in seq_along(maf_grid)) {
      wide[i, paste0("MDE_beta_MAF", maf_label(maf_grid[j]))] <- beta[j]
    }
    if (!all(is.na(beta_raw))) {
      for (j in seq_along(maf_grid)) {
        wide[i, paste0("MDE_beta_raw_MAF", maf_label(maf_grid[j]))] <- beta_raw[j]
      }
    }
    if (hazard) {
      for (j in seq_along(maf_grid)) {
        wide[i, paste0("MDE_HR_MAF", maf_label(maf_grid[j]))] <- hr[j]
      }
    }
  }

  long <- do.call(rbind, long)

  wide_path <- paste0(out_prefix, "_power_table.tsv")
  long_path <- paste0(out_prefix, "_power_long.tsv")
  write.table(round_numeric(wide), wide_path, sep = "\t", row.names = FALSE,
              quote = FALSE, na = "NA")
  write.table(round_numeric(long), long_path, sep = "\t", row.names = FALSE,
              quote = FALSE, na = "NA")
  invisible(list(wide = wide_path, long = long_path))
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
