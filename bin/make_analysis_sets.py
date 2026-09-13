#!/usr/bin/env python3
"""
Create analysis sets by extracting study arms and filtering samples.
Removes outliers and kinship-related samples from covariate data.
Used by MAKEANALYSISSETS process.
"""
import pandas as pd
import sys
import time
import re

from sample_ids import normalize_sample_ids


def normalize_to_iid_only(df):
    """
    Normalize table to IID-based format by removing FID/#FID and ensuring IID column name.
    """
    return normalize_sample_ids(df, source='covariates file')

def validate_column_names(df):
    """
    Validate column names for problematic characters.
    Returns error message if invalid characters found, None otherwise.
    """
    problematic_chars = r'[ +\-*/()[\]{}^~!@$%&=|\\<>?:;"\']'
    invalid_cols = []
    
    for col in df.columns:
        if re.search(problematic_chars, col):
            invalid_cols.append(col)
    
    if invalid_cols:
        error_msg = f"ERROR: Column names contain problematic characters (spaces, +, -, *, /, etc.):\n"
        for col in invalid_cols:
            error_msg += f"  - '{col}'\n"
        error_msg += "\nPlease use only alphanumeric characters and underscores in column names."
        return error_msg
    
    return None

def validate_cell_values(df):
    """
    Validate cell values for problematic characters (except #).
    Returns error message if invalid characters found, None otherwise.
    """
    # Allow # but not other special characters
    problematic_chars = r'[ +\*/()[\]{}^~!@$%&=|\\<>?:;"\']'
    invalid_entries = []
    
    for col in df.columns:
        # Only check string columns
        if df[col].dtype == 'object':
            for idx, value in df[col].items():
                if pd.notna(value) and isinstance(value, str):
                    if re.search(problematic_chars, value):
                        invalid_entries.append((col, idx, value))
                        if len(invalid_entries) >= 10:  # Limit to first 10 examples
                            break
        if len(invalid_entries) >= 10:
            break
    
    if invalid_entries:
        error_msg = f"ERROR: Cell values contain problematic characters (spaces, +, -, *, /, etc.):\n"
        error_msg += "First examples (column, row index, value):\n"
        for col, idx, value in invalid_entries[:10]:
            error_msg += f"  - Column '{col}', Row {idx}: '{value}'\n"
        if len(invalid_entries) > 10:
            error_msg += f"  ... and {len(invalid_entries) - 10} more\n"
        error_msg += "\nPlease use only alphanumeric characters, underscores, dots, and # in cell values."
        return error_msg
    
    return None

def main():
    if len(sys.argv) < 6:
        print("Usage: make_analysis_sets.py <samplelist.h5> <covariates.csv/tsv> <ancestry> <study_arm_col> <kinship_threshold>")
        sys.exit(1)
    
    samplelist_file = sys.argv[1]
    covariates_file = sys.argv[2]
    ancestry = sys.argv[3]
    study_arm_colname = sys.argv[4]
    kinship_threshold = float(sys.argv[5])
    
    # Read data files
    ancestry_df = pd.read_hdf(samplelist_file, key="ancestry_keep")
    outlier_df = pd.read_hdf(samplelist_file, key="outliers")
    kin_df = pd.read_hdf(samplelist_file, key="kin")
    
    # Convert IID columns to strings for consistent comparison
    # HDF5 files may store IIDs as integers, but we need them as strings to match covariates
    for df in [ancestry_df, outlier_df]:
        if 'IID' in df.columns:
            df['IID'] = df['IID'].astype(str)
        if '#IID' in df.columns:
            df['#IID'] = df['#IID'].astype(str)
    
    # For kinship dataframe, convert both IID columns
    if not kin_df.empty:
        if '#IID1' in kin_df.columns:
            kin_df['#IID1'] = kin_df['#IID1'].astype(str)
        if 'IID2' in kin_df.columns:
            kin_df['IID2'] = kin_df['IID2'].astype(str)
    
    # Auto-detect delimiter (CSV or TSV)
    with open(covariates_file, 'r') as f:
        first_line = f.readline()
        delimiter = '\t' if '\t' in first_line else ','
    
    # Read IDs as strings and normalize to IID-only format for downstream PLINK2 compatibility.
    data_df = pd.read_csv(
        covariates_file,
        sep=delimiter,
        engine='c',
        dtype={'#FID': str, 'FID': str, '#IID': str, 'IID': str}
    )
    data_df = normalize_to_iid_only(data_df)
    
    # Validate column names
    validation_error = validate_column_names(data_df)
    if validation_error:
        print(validation_error, file=sys.stderr)
        sys.exit(1)
    
    # Validate cell values
    validation_error = validate_cell_values(data_df)
    if validation_error:
        print(validation_error, file=sys.stderr)
        sys.exit(1)
    
    # Apply filtering steps before splitting by study arm
    print('=' * 80)
    print('FILTERING ANALYTICAL SET')
    print('=' * 80)
    
    initial_count = len(data_df)
    # error if covariates files has non-unique IIDs
    if data_df.IID.duplicated().any():
        dup_iids = data_df.IID[data_df.IID.duplicated()].unique()
        error_msg = f"ERROR: Covariates file contains non-unique IIDs:\n"
        for iid in dup_iids:
            error_msg += f"  - {iid}\n"
        print(error_msg, file=sys.stderr)
        sys.exit(1)
    print(f'Initial participants count in {covariates_file}: {initial_count}')
    
    # Step 1: Filter by ancestry
    analytical_set = data_df[data_df.IID.isin(ancestry_df.IID)].copy(deep=True)
    ancestry_filtered_count = len(analytical_set)
    print(f'After ancestry filtering: {ancestry_filtered_count} (removed {initial_count - ancestry_filtered_count})')
    
    # Step 2: Remove outliers
    if not outlier_df.empty:
        outlier_count = analytical_set.IID.isin(outlier_df.IID).sum()
        analytical_set = analytical_set[~analytical_set.IID.isin(outlier_df.IID)].copy(deep=True)
        print(f'After outlier removal: {len(analytical_set)} (removed {outlier_count})')
    else:
        print(f'No outliers to remove')
    
    # Step 3: Remove kinship-related samples
    kin_df = kin_df[kin_df.KINSHIP >= kinship_threshold]
    kin_pairs = kin_df[(kin_df['#IID1'].isin(analytical_set.IID)) & 
                       (kin_df.IID2.isin(analytical_set.IID))].copy()
    kinship_removed_count = kin_pairs.shape[0]
    analytical_set = analytical_set[~analytical_set.IID.isin(kin_pairs.IID2)].copy()
    print(f'After kinship removal: {len(analytical_set)} (removed {kinship_removed_count})')
    
    # Write out the full analytical set (all study arms combined)
    analytical_set_filename = f"{ancestry}_all.tsv"
    analytical_set.to_csv(analytical_set_filename, sep="\t", index=False)
    print(f'\nFull analytical set saved: {analytical_set_filename}')
    print(f'Final analytical set: {len(analytical_set)} samples')
    
    # Extract all unique study arms
    study_arms = analytical_set[study_arm_colname].unique().tolist()
    print(f'\nFound {len(study_arms)} study arms in analytical set: {study_arms}')
    print('=' * 80)
    
    # Process each study arm from the pre-filtered analytical set
    for study_arm in study_arms:
        print(f'\n---- {study_arm} ----')
        
        study_arm_samples = analytical_set[analytical_set[study_arm_colname] == study_arm].copy(deep=True)
        sample_count = len(study_arm_samples)
        print(f'Samples in {study_arm}: {sample_count}')
        
        # Enforce minimum sample requirement
        if sample_count <= 50:
            error_msg = f"\nERROR: Insufficient samples in study arm '{study_arm}'\n"
            error_msg += f"Found {sample_count} samples, but minimum required is 51.\n"
            error_msg += f"This study arm/type does not meet the minimum sample requirement for analysis.\n"
            print(error_msg, file=sys.stderr)
            sys.exit(1)
        
        study_arm_samples.to_csv(f"{ancestry}_{study_arm}_filtered.tsv", sep="\t", index=False)
    
    time.sleep(5)

if __name__ == "__main__":
    main()
