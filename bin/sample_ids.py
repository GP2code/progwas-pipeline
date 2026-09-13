"""
Shared sample-identifier handling for user-supplied phenotype/covariate tables.

The pipeline keys samples on IID alone. Genotype files may carry an FID
(PLINK1 .fam always does, a .psam may), but that is handled separately: the
genetic-QC steps run bin/normalize_psam_iid_only.sh, and bin/make_keep_iid.py
reconstructs the FID from the .psam when a plink2 --keep file is needed.

Phenotype and covariate files are different. An FID column there is accepted but
carries no meaning for the pipeline, so it is dropped on read. That keeps the
two accepted input shapes -- "#FID IID ..." and "IID ..." -- behaviourally
identical, which is why the smoke-test matrix does not need to cover both shapes
once per model.

Accepted header spellings for the identifier: IID, #IID.
Accepted (and ignored) family-identifier spellings: FID, #FID.
"""

IID_ALIASES = ('IID', '#IID')
FID_ALIASES = ('FID', '#FID')


def normalize_sample_ids(df, source='table', copy=True):
    """Return `df` keyed on a plain string `IID` column, with any FID dropped.

    Raises ValueError if no identifier column is present.
    """
    out = df.copy(deep=True) if copy else df

    # '#IID' is the plink2 spelling; fold it into a plain 'IID'.
    if 'IID' not in out.columns:
        for alias in IID_ALIASES:
            if alias in out.columns:
                out.rename(columns={alias: 'IID'}, inplace=True)
                break

    if 'IID' not in out.columns:
        raise ValueError(
            f"'IID' column not found in {source}. "
            f"Found columns: {list(out.columns)}. "
            "Phenotype and covariate files must have an 'IID' (or '#IID') column."
        )

    dropped = [c for c in FID_ALIASES if c in out.columns]
    if dropped:
        out.drop(columns=dropped, inplace=True)
        print(f"Ignoring family-ID column(s) {dropped} in {source}; "
              "samples are matched on IID only.")

    out['IID'] = out['IID'].astype(str).str.strip()
    return out
