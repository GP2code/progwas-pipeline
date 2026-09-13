#!/usr/bin/env python3
"""
Preprocess phenotype and covariate data for PLINK GLM analysis.
Handles standardization, one-hot encoding, and phenotype format conversion.
"""

import pandas as pd
import numpy as np
import sys
import argparse

from sample_ids import normalize_sample_ids


def main():
    parser = argparse.ArgumentParser(description='Preprocess phenotype and covariate data')
    parser.add_argument('--samplelist', required=True, help='Sample list file with PCA data')
    parser.add_argument('--phenofile', required=True, help='Phenotype file')
    parser.add_argument('--outfile', required=True, help='Output file prefix')
    parser.add_argument('--pheno-name', default='', help='Space-separated phenotype names')
    parser.add_argument('--covar-numeric', default='', help='Space-separated numeric covariate names')
    parser.add_argument('--covar-categorical', default='', help='Space-separated categorical covariate names')
    parser.add_argument('--covar-interact', default='', help='Interaction covariate (not standardized)')
    
    args = parser.parse_args()
    
    # Create log file and redirect stdout
    log_filename = f"{args.outfile}_preprocessing.log"
    log_file = open(log_filename, "w", buffering=1)
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    
    try:
        sys.stdout = log_file
        sys.stderr = log_file
        
        print("=" * 80)
        print(f"EXPORT_PLINK Preprocessing Log")
        print(f"Study arm: {args.outfile}")
        print(f"Timestamp: {pd.Timestamp.now()}")
        print("=" * 80)
        print()

        all_phenos = args.pheno_name.split() if args.pheno_name else []
        covar_num = args.covar_numeric.split() if args.covar_numeric else []
        covar_cat = args.covar_categorical.split() if args.covar_categorical else []
        interact_covar = args.covar_interact.strip()
        
        with open(args.phenofile, 'r') as f:
            pheno_delim = '\t' if '\t' in f.readline() else ','
        d_pheno = pd.read_csv(args.phenofile, sep=pheno_delim, engine='c')
        d_sample = pd.read_csv(args.samplelist, sep="\t", engine='c')

        # Key on IID and discard any family-ID column; see bin/sample_ids.py.
        d_pheno = normalize_sample_ids(d_pheno, source=args.phenofile)
        d_sample = normalize_sample_ids(d_sample, source=args.samplelist)

        # Drop columns from d_sample that already exist in d_pheno (d_pheno takes precedence,
        # as it may contain time-varying variables that differ from the static covariate file)
        overlap_cols = [c for c in d_pheno.columns if c in d_sample.columns and c != 'IID']
        if overlap_cols:
            print(f"Dropping overlapping columns from sample list (kept from phenotype file): {overlap_cols}")
            d_sample = d_sample.drop(columns=overlap_cols)

        d_result = pd.merge(d_pheno, d_sample, on='IID', how='inner')
        if d_result.shape[0] == 0:
            print("WARNING: No samples found after merging phenotype and sample files")
            print("=" * 80)
            print("Preprocessing completed with warnings")
            print("=" * 80)
            log_file.close()
            sys.stdout = original_stdout
            sys.stderr = original_stderr
            sys.exit(0)
        
        # Only include covariates explicitly specified by user
        d_set = d_result.loc[:, ["IID"] + all_phenos + covar_num].copy()
        print(f"Including numeric covariates: {covar_num}")
        
        # One-hot encode categorical covariates
        categorical_dummies = []
        if covar_cat:
            for cat_col in covar_cat:
                if cat_col in d_result.columns:
                    # Use most frequent level as reference to stabilize design matrix.
                    cat_series = d_result[cat_col].astype('string')
                    level_counts = cat_series.dropna().value_counts()
                    if level_counts.empty:
                        print(f"Skipping '{cat_col}' (no non-missing values)")
                        continue
                    ref_level = level_counts.idxmax()
                    other_levels = sorted([x for x in level_counts.index.tolist() if x != ref_level])
                    ordered_levels = [ref_level] + other_levels
                    cat_series = pd.Categorical(cat_series, categories=ordered_levels, ordered=True)

                    # Create dummy variables and drop the explicit reference level.
                    dummies = pd.get_dummies(cat_series, prefix=cat_col, drop_first=True, dtype=int)
                    d_set = pd.concat([d_set, dummies], axis=1)
                    categorical_dummies.extend(list(dummies.columns))
                    print(
                        f"One-hot encoded '{cat_col}' with reference='{ref_level}': {list(dummies.columns)}"
                    )
        
        # Standardize numeric covariates (except interaction covariate and categorical dummies)
        print("\n--- Standardizing numeric covariates ---")
        for col in covar_num:
            if col in d_set.columns:
                # The interaction covariate is CENTERED but not scaled.
                #
                # Not scaling keeps the interaction coefficient in the covariate's
                # original units (e.g. per year of age), which is the point of
                # excluding it here. Centering is still required: plink2 --glm
                # interaction builds an explicit SNP x covar column, and if the
                # covariate has a large non-zero mean that column is nearly a
                # multiple of the SNP column. The resulting collinearity trips
                # plink2's variance-inflation gate (default --vif 50) and every
                # affected variant is dropped with ERRCODE=VIF_TOO_HIGH rather
                # than analysed -- on the bundled test data, an uncentered
                # age_at_baseline (mean 72.9) lost 1525 of 1558 variants.
                #
                # Centering leaves the interaction estimate unchanged and merely
                # redefines the SNP main effect as the effect at the mean
                # covariate value, which is the usual parameterisation anyway.
                if interact_covar and col == interact_covar:
                    if pd.api.types.is_numeric_dtype(d_set[col]):
                        covar_mean = d_set[col].mean()
                        d_set[col] = d_set[col] - covar_mean
                        print(f"Centered interaction covariate '{col}' (mean={covar_mean:.3f}); "
                              "kept on original scale for interaction interpretation")
                    else:
                        print(f"Keeping '{col}' on original scale for interaction interpretation")
                    continue

                # Check if numeric
                if pd.api.types.is_numeric_dtype(d_set[col]):
                    unique_vals = d_set[col].dropna().unique()
                    # Standardize if not a binary 0/1 variable
                    if not (set(unique_vals).issubset({0, 1}) and len(unique_vals) == 2):
                        # Standardize: mean=0, variance=1 (using unbiased variance estimator)
                        covar_mean = d_set[col].mean()
                        covar_std = d_set[col].std(ddof=1)
                        if covar_std > 0:
                            d_set[col] = (d_set[col] - covar_mean) / covar_std
                            print(f"Standardized covariate: {col} (mean={covar_mean:.3f}, sd={covar_std:.3f})")
                        else:
                            print(f"Skipping standardization for {col} (constant variable)")

        for col in covar_num:
            if col in d_set.columns and d_set[col].isna().any():
                n_missing = d_set[col].isna().sum()
                d_set[col] = d_set[col].astype(object).where(d_set[col].notna(), 'NA')
                print(f"Numeric covariate '{col}': filled {n_missing} missing values with 'NA'")
        
        # Skip standardization for categorical dummy variables (they're already 0/1)
        for dummy_col in categorical_dummies:
            if dummy_col in d_set.columns:
                print(f"Skipping standardization for categorical dummy variable: {dummy_col}")
        
        # Convert binary phenotypes (0/1/missing) to PLINK format (1/2/-9)
        for pheno_col in all_phenos:
            if pheno_col in d_set.columns:
                unique_vals = d_set[pheno_col].dropna().unique()
                # Check if it's binary 0/1 coded (not already 1/2/-9)
                if set(unique_vals).issubset({0, 1}) and len(unique_vals) > 0:
                    print(f"Converting binary phenotype '{pheno_col}' from 0/1 to PLINK format 1/2/-9")
                    # 0 -> 1 (control), 1 -> 2 (case), NaN -> -9 (missing)
                    d_set[pheno_col] = d_set[pheno_col].map({0: 1, 1: 2}).fillna(-9).astype(int)
                elif set(unique_vals).issubset({1, 2, -9}):
                    # Already in PLINK format, ensure missing is -9
                    d_set[pheno_col] = d_set[pheno_col].fillna(-9).astype(int)
                    print(f"Phenotype '{pheno_col}' already in PLINK format 1/2/-9")
                elif d_set[pheno_col].isna().any():
                    n_missing = d_set[pheno_col].isna().sum()
                    d_set[pheno_col] = d_set[pheno_col].astype(object).where(d_set[pheno_col].notna(), 'NA')
                    print(f"Quantitative phenotype '{pheno_col}': filled {n_missing} missing values with 'NA'")
        
        # Extract all covariate column names (numeric + categorical dummies)
        all_covar_cols = [col for col in d_set.columns if col not in ["IID"] + all_phenos]
        
        # Reorder covariates if interaction term is specified
        if interact_covar:
            if interact_covar in all_covar_cols:
                # Put interaction covariate first
                all_covar_cols.remove(interact_covar)
                all_covar_cols = [interact_covar] + all_covar_cols
                print(f"\nInteraction analysis: '{interact_covar}' as first covariate")
            else:
                print(f"\nWARNING: Interaction covariate '{interact_covar}' not found in data")
        
        # Write covariate metadata
        covar_names_str = ','.join(all_covar_cols)
        with open(f"{args.outfile}_covar_names.txt", 'w') as f:
            f.write(covar_names_str)
        
        with open(f"{args.outfile}_n_covar.txt", 'w') as f:
            f.write(str(len(all_covar_cols)))
        
        print(f"\nTotal covariates: {len(all_covar_cols)}")
        print(f"Covariate order: {covar_names_str}")
        
        # Save single file with all phenotypes (plink2 --glm handles missing phenotypes automatically)
        output_path = f"{args.outfile}_filtered.pca.pheno.tsv"
        d_set.to_csv(output_path, sep="\t", index=False)

        with open(output_path) as f:
            header_fields = f.readline().rstrip('\n').split('\t')
            expected_n = len(header_fields)
            for lineno, line in enumerate(f, start=2):
                n_fields = len(line.rstrip('\n').split('\t'))
                if n_fields != expected_n:
                    raise ValueError(
                        f"{output_path} line {lineno} has {n_fields} tab-separated fields, "
                        f"expected {expected_n} (header: {header_fields}). Offending line (truncated): {line[:200]!r}"
                    )
        print(f"\nCreated {args.outfile}_filtered.pca.pheno.tsv with {d_set.shape[0]} samples")
        print(f"Phenotypes: {', '.join(all_phenos)}")
        for pheno_col in all_phenos:
            if pheno_col in d_set.columns:
                # Count valid values: exclude -9 for binary phenotypes in PLINK format
                unique_vals = d_set[pheno_col].dropna().unique()
                if set(unique_vals).issubset({1, 2, -9}):
                    # Binary phenotype in PLINK format: count only 1 and 2
                    n_valid = ((d_set[pheno_col] == 1) | (d_set[pheno_col] == 2)).sum()
                else:
                    # Quantitative phenotype: count non-missing while excluding explicit NA tokens
                    n_valid = (d_set[pheno_col].notna() & (d_set[pheno_col] != 'NA')).sum()
                print(f"  {pheno_col}: {n_valid} samples with valid values")
        
        print()
        print("=" * 80)
        print("Preprocessing completed successfully")
        print("=" * 80)
        
        # Close log file and restore stdout/stderr
        log_file.close()
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        
    except Exception as e:
        # Ensure we close log and restore stdout/stderr even on error
        log_file.close()
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        print(f"Error during preprocessing: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        raise


if __name__ == '__main__':
    main()
