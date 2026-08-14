#!/usr/bin/env python3
"""
Manhattan and QQ plot generator for GWAS results

Usage:
    manhattan.py --input <results.tsv> --model <glm|cph|lmm_gallop> [--suffix <str>]

Arguments:
    --input: Path to GWAS results TSV file
    --model: Analysis model type (glm, cph, or lmm_gallop)
    --suffix: Optional suffix for cohort naming (default: "null")
"""

import argparse
import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
from qmplot import manhattanplot, qqplot


def plot_summary_stats(data, cohort, outcome, model, p_col='P'):
    """Generate Manhattan and QQ plots for GWAS summary statistics.
    
    Args:
        data: pandas DataFrame with GWAS results
        cohort: Cohort name for plot titles and filenames
        outcome: Phenotype/outcome name
        model: Model type (glm, cph, or lmm_gallop)
    """
    xtick = set(['chr' + i for i in list(map(str, range(1, 14))) + ['15', '17', '19', '22']])

    if model == "lmm_gallop":
        # GALLOP longitudinal model - plot both intercept and slope
        
        # Intercept Manhattan plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        manhattanplot(data=data,
                      title=f"Manhattan Intercept {cohort} {outcome}",
                      pv="P", ax=ax, 
                      xtick_label_set=xtick)
        plt.savefig(f"{cohort}_{outcome}_manhattan_intercept.{model}.png", dpi=300)
        
        # Intercept QQ plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        qqplot(data=data["P"],
               marker="o",
               title=f"QQ Intercept {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               ax=ax)
        plt.savefig(f"{cohort}_{outcome}_qq_intercept.{model}.png", dpi=300)
        
        # Slope Manhattan plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        manhattanplot(data=data,
                      title=f"Manhattan Slope {cohort} {outcome}",
                      pv="P_INT", ax=ax, 
                      xtick_label_set=xtick)
        plt.savefig(f"{cohort}_{outcome}_manhattan_slope.{model}.png", dpi=300)
        
        # Slope QQ plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        qqplot(data=data["P_INT"],
               marker="o",
               title=f"QQ Slope {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               ax=ax)
        plt.savefig(f"{cohort}_{outcome}_qq_slope.{model}.png", dpi=300)
    
    elif model in ["glm", "cph"]:
        # Standard GLM or Cox PH model - single P-value
        
        # Manhattan plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        manhattanplot(data=data,
                      title=f"Manhattan {model} {cohort} {outcome}",
                      pv=p_col, ax=ax,
                      xtick_label_set=xtick)
        plt.savefig(f"{cohort}_{outcome}_manhattan.{model}.png", dpi=300)
        
        # QQ plot
        f, ax = plt.subplots(figsize=(15, 7), facecolor="w", edgecolor="k")
        qqplot(data=data[p_col],
               marker="o",
               title=f"QQ {model} {cohort} {outcome}",
               xlabel=r"Expected -log(P)",
               ylabel=r"Observed -log(P)",
               ax=ax)
        plt.savefig(f"{cohort}_{outcome}_qq.{model}.png", dpi=300)
    
    else:
        raise ValueError(f"Model '{model}' not recognized. Use 'glm', 'cph', or 'lmm_gallop'.")
    
    print(f"✅ Results plotting success for {cohort}_{outcome} ({model})")


def main():
    """Main entry point for Manhattan plot generation."""
    parser = argparse.ArgumentParser(
        description="Generate Manhattan and QQ plots from GWAS results",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input", required=True, help="Input TSV file with GWAS results")
    parser.add_argument("--model", required=True, choices=["glm", "cph", "lmm_gallop"],
                       help="Analysis model type")
    parser.add_argument("--p-col", default="P",
                       help="P-value column for glm/cph plots (default: P)")
    parser.add_argument("--suffix", default="null", help="Cohort suffix (default: null)")
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"❌ Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Extract metadata from filename
    # Expected format: {ancestry}_{cohort}_{phenotype}_allresults.tsv
    gwas_name = os.path.splitext(os.path.basename(args.input))[0]
    
    try:
        cohort = gwas_name.split('_')[1]
        find_phenoname = cohort + r'_(.*?)_allresults'
        result = re.search(find_phenoname, gwas_name)
        if not result:
            raise ValueError(f"Could not extract phenotype name from filename: {gwas_name}")
        pheno = result.group(1)
    except (IndexError, ValueError) as e:
        print(f"❌ Error parsing filename: {e}", file=sys.stderr)
        print(f"Expected format: <ancestry>_<cohort>_<phenotype>_allresults.tsv", file=sys.stderr)
        sys.exit(1)
    
    # Add suffix to cohort if provided
    cohort_suffix = f"{cohort}.{args.suffix}" if args.suffix and args.suffix != "null" else cohort
    
    # Load and prepare data
    print(f"📊 Loading data from: {args.input}")
    try:
        df = pd.read_csv(args.input, sep="\t")
    except Exception as e:
        print(f"❌ Error reading file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"   Loaded {len(df)} variants")
    
    # Drop rows missing the columns required for this model's plots
    required_cols = ['#CHROM', 'POS']
    if args.model in ['glm', 'cph']:
        required_cols.append('P')
    elif args.model == 'lmm_gallop':
        required_cols.extend(['P', 'P_INT'])
    df = df.dropna(subset=required_cols)
    print(f"   After removing missing values: {len(df)} variants")
    
    # Add chromosome ordering for proper sorting
    df['chr_order'] = df['#CHROM'].str.replace('chr', '')
    df['chr_order'] = df['chr_order'].astype(int)
    df = df.sort_values(by=['chr_order', 'POS'])
    
    # Generate plots
    print(f"🎨 Generating plots for cohort={cohort_suffix}, phenotype={pheno}, model={args.model}, p_col={args.p_col}")
    plot_summary_stats(data=df, cohort=cohort_suffix, outcome=pheno, model=args.model, p_col=args.p_col)


if __name__ == "__main__":
    main()
