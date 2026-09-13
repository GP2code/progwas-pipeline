#!/usr/bin/env python3
"""
Generate Table 1 (descriptive statistics) and Kaplan-Meier curves for GWAS analyses.
Reads configuration from YAML file and produces:
- Table 1 CSV with descriptive statistics by study arm
- KM survival curves (if survival_flag is true)
"""

import argparse
import pandas as pd
import numpy as np
import sys
from tableone import TableOne
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for server environments
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

from sample_ids import normalize_sample_ids

def main():
    parser = argparse.ArgumentParser(
        description='Generate Table 1 and Kaplan-Meier plots for GWAS analyses'
    )
    parser.add_argument('covarfile',      help='Covariate file (TSV)')
    parser.add_argument('phenofile',      help='Phenotype file (TSV)')
    parser.add_argument('analytical_set', help='Filtered analysis set from MAKEANALYSISSETS (TSV)')
    parser.add_argument('--pheno-name',        required=True,         dest='pheno_name')
    parser.add_argument('--study-arm-col',     default='study_arm',   dest='study_arm_col')
    parser.add_argument('--covar-numeric',     nargs='*', default=[], dest='covar_numeric')
    parser.add_argument('--covar-categorical', nargs='*', default=[], dest='covar_categorical')
    parser.add_argument('--time-col',          default='study_days',  dest='time_col')
    parser.add_argument('--longitudinal',      action='store_true')
    parser.add_argument('--survival',          action='store_true')
    parser.add_argument('--ancestry',          default='EUR')
    parser.add_argument('--analysis-name',     default='TEST',        dest='analysis_name')
    args = parser.parse_args()

    covarfile         = args.covarfile
    phenofile         = args.phenofile
    analytical_set    = args.analytical_set
    pheno_name        = args.pheno_name
    study_arm_col     = args.study_arm_col
    covar_numeric     = args.covar_numeric
    covar_categorical = args.covar_categorical
    time_col          = args.time_col
    longitudinal_flag = args.longitudinal
    survival_flag     = args.survival
    ancestry          = args.ancestry
    analysis_name     = args.analysis_name

    print(f"Analysis: {analysis_name}")
    print(f"Ancestry: {ancestry}")
    print()

    # Read covariate and phenotype files
    print(f"Reading covariate file: {covarfile}")
    with open(covarfile, 'r') as f:
        cov_delim = '\t' if '\t' in f.readline() else ','
    cov = pd.read_csv(covarfile, sep=cov_delim)
    cov = normalize_sample_ids(cov, source=covarfile)
    print(f"Covariates shape: {cov.shape}")

    print(f"Reading phenotype file: {phenofile}")
    with open(phenofile, 'r') as f:
        pheno_delim = '\t' if '\t' in f.readline() else ','
    pheno = pd.read_csv(phenofile, sep=pheno_delim)
    pheno = normalize_sample_ids(pheno, source=phenofile)
    print(f"Phenotype shape: {pheno.shape}")
    print(f"Phenotype column: {pheno_name}")
    print()

    print(f"Reading filtered analysis set: {analytical_set}")
    filtered_covar = pd.read_csv(analytical_set, sep='\t')
    print(f"Samples after filtering: {filtered_covar.shape[0]}")

    # Merge covariates and phenotypes for Table 1
    print("=" * 80)
    print("GENERATING TABLE 1 - Descriptive Statistics by Study Arm")
    print("=" * 80)

    # Merge on IID
    merged = pd.merge(cov, pheno, on='IID', how='inner', suffixes=('_cov', '_pheno'))
    merged = merged[merged['IID'].isin(filtered_covar['IID'])]
    print(f"Merged data shape: {merged.shape}")

    # Handle longitudinal data if needed
    if longitudinal_flag:
        print("\nLongitudinal analysis detected - processing repeated measures...")
        
        # Identify time column (handle suffix conflicts)
        time_col_actual = time_col
        if time_col not in merged.columns:
            if f"{time_col}_cov" in merged.columns:
                time_col_actual = f"{time_col}_cov"
            elif f"{time_col}_pheno" in merged.columns:
                time_col_actual = f"{time_col}_pheno"
        
        # Count observations per IID
        obs_counts = merged.groupby('IID').size().reset_index(name='N_obs')
        
        # Get last follow-up time per IID
        last_obs = merged.groupby('IID')[time_col_actual].max().reset_index()
        last_obs.columns = ['IID', 'last_obs_time']
        
        # Get first observation per IID for Table 1
        merged_first = merged.sort_values(['IID', time_col_actual]).groupby('IID').first().reset_index()
        
        # Merge N_obs and last_obs_time
        merged_first = merged_first.merge(obs_counts, on='IID', how='left')
        merged_first = merged_first.merge(last_obs, on='IID', how='left')
        
        print(f"Total observations: {len(merged)}")
        print(f"Unique subjects: {len(merged_first)}")
        print(f"Median observations per subject: {obs_counts['N_obs'].median():.1f}")
        print(f"Median last follow-up: {last_obs['last_obs_time'].median():.1f}")
        
        # Use first observation for Table 1
        merged = merged_first
        print(f"Using first observation per subject for Table 1: {merged.shape}")
        print()

    print(f"Final merged data shape: {merged.shape}")

    # Identify study_arm column (handle suffix conflicts)
    if study_arm_col in merged.columns:
        arm_col = study_arm_col
    elif f"{study_arm_col}_cov" in merged.columns:
        arm_col = f"{study_arm_col}_cov"
    elif f"{study_arm_col}_pheno" in merged.columns:
        arm_col = f"{study_arm_col}_pheno"
    else:
        print(f"Warning: study_arm column '{study_arm_col}' not found. Available columns: {list(merged.columns)}")
        arm_col = None

    if arm_col:
        study_arms = sorted(merged[arm_col].unique())
        print(f"Study arm column: {arm_col}")
        print(f"Study arms: {study_arms}")
        print(f"Number of unique study arms: {len(study_arms)}")
        print()

        # Prepare columns for Table 1
        columns_for_table = []
        categorical_cols = []
        nonnormal_cols = []  # Columns to display as median [Q1, Q3]
        
        # Auto-detect if phenotype is categorical (< 5 unique values)
        pheno_col_actual = pheno_name
        if pheno_name not in merged.columns:
            if f"{pheno_name}_cov" in merged.columns:
                pheno_col_actual = f"{pheno_name}_cov"
            elif f"{pheno_name}_pheno" in merged.columns:
                pheno_col_actual = f"{pheno_name}_pheno"
        
        if pheno_col_actual in merged.columns:
            n_unique = merged[pheno_col_actual].nunique()
            print(f"Phenotype '{pheno_name}' has {n_unique} unique values")
            
            columns_for_table.append(pheno_col_actual)
            if n_unique < 5:
                categorical_cols.append(pheno_col_actual)
                print(f"  → Treating as categorical/binomial")
            else:
                print(f"  → Treating as continuous")
        
        # Add time column for survival/longitudinal analysis (use median [Q1, Q3])
        time_col_actual = time_col
        if time_col not in merged.columns:
            if f"{time_col}_cov" in merged.columns:
                time_col_actual = f"{time_col}_cov"
            elif f"{time_col}_pheno" in merged.columns:
                time_col_actual = f"{time_col}_pheno"
        
        if time_col_actual in merged.columns:
            columns_for_table.append(time_col_actual)
            nonnormal_cols.append(time_col_actual)  # Display as median [Q1, Q3]
        
        # Add longitudinal-specific columns
        if longitudinal_flag:
            if 'N_obs' in merged.columns:
                columns_for_table.append('N_obs')
                nonnormal_cols.append('N_obs')
            if 'last_obs_time' in merged.columns:
                columns_for_table.append('last_obs_time')
                nonnormal_cols.append('last_obs_time')
        
        # Add numeric covariates
        for col in covar_numeric:
            if col in merged.columns:
                columns_for_table.append(col)
            elif f"{col}_cov" in merged.columns:
                columns_for_table.append(f"{col}_cov")
            elif f"{col}_pheno" in merged.columns:
                columns_for_table.append(f"{col}_pheno")
        
        # Add categorical covariates
        for col in covar_categorical:
            if col in merged.columns:
                columns_for_table.append(col)
                categorical_cols.append(col)
            elif f"{col}_cov" in merged.columns:
                columns_for_table.append(f"{col}_cov")
                categorical_cols.append(f"{col}_cov")
            elif f"{col}_pheno" in merged.columns:
                columns_for_table.append(f"{col}_pheno")
                categorical_cols.append(f"{col}_pheno")
        
        # Remove duplicates while preserving order
        columns_for_table = list(dict.fromkeys(columns_for_table))
        
        # Filter to only existing columns
        columns_for_table = [col for col in columns_for_table if col in merged.columns]
        
        # Remove columns with no variation (would cause statistical test failures)
        no_variation = [col for col in columns_for_table if merged[col].nunique() <= 1]
        if no_variation:
            print(f"Dropping zero-variation columns: {no_variation}")
        columns_for_table = [col for col in columns_for_table if col not in no_variation]
        categorical_cols = [col for col in categorical_cols if col not in no_variation]
        nonnormal_cols = [col for col in nonnormal_cols if col not in no_variation]
        
        print(f"Variables for Table 1: {columns_for_table}")
        print(f"Categorical variables: {categorical_cols}")
        print(f"Non-normal variables (median [Q1, Q3]): {nonnormal_cols}")
        print()
        
        # Generate Table 1
        try:
            if len(study_arms) >= 2:
                # Multiple study arms - generate grouped table with comparisons
                table1 = TableOne(
                    merged,
                    columns=columns_for_table,
                    categorical=categorical_cols,
                    nonnormal=nonnormal_cols,
                    groupby=arm_col,
                    pval=True,
                    missing=True
                )
            else:
                # Single study arm - generate overall descriptive statistics only
                print("Note: Only one study arm found. Generating overall descriptive statistics.")
                table1 = TableOne(
                    merged,
                    columns=columns_for_table,
                    categorical=categorical_cols,
                    nonnormal=nonnormal_cols,
                    pval=False,
                    missing=True
                )
            
            print(table1)
            print()
            
            # Save to current directory (Nextflow publishDir handles routing)
            output_file = f"table1_{analysis_name}_{pheno_name}.csv"
            table1.to_csv(output_file)
            print(f"Table 1 saved to: {output_file}")
            print()
        except Exception as e:
            print(f"Error generating Table 1: {e}")
            import traceback
            traceback.print_exc()
            print()
    else:
        print("Skipping Table 1 generation due to missing study_arm column.")
        print()

    # ==================================================================================
    # KAPLAN-MEIER SURVIVAL CURVES (if survival_flag is true)
    # ==================================================================================
    if survival_flag and arm_col:
        print("=" * 80)
        print("GENERATING KAPLAN-MEIER SURVIVAL CURVES")
        print("=" * 80)
        
        try:
            # Reload full phenotype data (before first-obs reduction)
            merged_full = pd.merge(cov, pheno, on='IID', how='inner', suffixes=('_cov', '_pheno'))
            merged_full = merged_full[merged_full['IID'].isin(filtered_covar['IID'])]
            
            # Identify time and event columns
            time_col_actual = time_col
            if time_col not in merged_full.columns:
                if f"{time_col}_cov" in merged_full.columns:
                    time_col_actual = f"{time_col}_cov"
                elif f"{time_col}_pheno" in merged_full.columns:
                    time_col_actual = f"{time_col}_pheno"
            
            pheno_col_actual = pheno_name
            if pheno_name not in merged_full.columns:
                if f"{pheno_name}_cov" in merged_full.columns:
                    pheno_col_actual = f"{pheno_name}_cov"
                elif f"{pheno_name}_pheno" in merged_full.columns:
                    pheno_col_actual = f"{pheno_name}_pheno"
            
            arm_col_actual = arm_col
            if arm_col not in merged_full.columns:
                if f"{arm_col}_cov" in merged_full.columns:
                    arm_col_actual = f"{arm_col}_cov"
                elif f"{arm_col}_pheno" in merged_full.columns:
                    arm_col_actual = f"{arm_col}_pheno"
            
            # Check if we have the required columns
            if time_col_actual in merged_full.columns and pheno_col_actual in merged_full.columns:
                # Remove missing values
                km_data = merged_full[[arm_col_actual, time_col_actual, pheno_col_actual]].dropna()
                
                print(f"Time column: {time_col_actual}")
                print(f"Event column: {pheno_col_actual}")
                print(f"Study arm column: {arm_col_actual}")
                print(f"Sample size for KM analysis: {len(km_data)}")
                print()
                
                # Create figure
                plt.figure(figsize=(10, 6))
                
                # Fit and plot KM curve for each study arm
                kmf = KaplanMeierFitter()
                study_arms_km = sorted(km_data[arm_col_actual].unique())
                
                for arm in study_arms_km:
                    arm_data = km_data[km_data[arm_col_actual] == arm]
                    n_events = int(arm_data[pheno_col_actual].sum())
                    n_total = len(arm_data)
                    
                    kmf.fit(
                        durations=arm_data[time_col_actual],
                        event_observed=arm_data[pheno_col_actual],
                        label=f'Arm {arm} (n={n_total}, events={n_events})'
                    )
                    kmf.plot_survival_function(ci_show=True)
                    
                    print(f"Study arm {arm}: n={n_total}, events={n_events} ({100*n_events/n_total:.1f}%)")
                
                plt.xlabel(f'Time ({time_col})')
                plt.ylabel('Survival Probability')
                plt.title(f'Kaplan-Meier Survival Curves by Study Arm\n{analysis_name}')
                plt.legend(loc='best')
                plt.grid(True, alpha=0.3)
                
                # Save to current directory (Nextflow publishDir handles routing)
                km_plot_file = f"{ancestry}_{pheno_name}_km_plot.png"
                plt.savefig(km_plot_file, dpi=300, bbox_inches='tight')
                print(f"\nKaplan-Meier plot saved to: {km_plot_file}")
                plt.close()
                
            else:
                print(f"Warning: Could not find required columns for KM analysis")
                print(f"  Time column '{time_col}' found as: {time_col_actual if time_col_actual in merged_full.columns else 'NOT FOUND'}")
                print(f"  Event column '{pheno_name}' found as: {pheno_col_actual if pheno_col_actual in merged_full.columns else 'NOT FOUND'}")
        
        except Exception as e:
            print(f"Error generating Kaplan-Meier curves: {e}")
            import traceback
            traceback.print_exc()
        
        print()

    print("=" * 80)
    print("TABLE 1 GENERATION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
