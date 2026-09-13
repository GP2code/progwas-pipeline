# Pipeline smoke test

`tests/smoke_test.sh` runs the main pipeline (`main.nf`) end to end over the
bundled example data (chr20–22, ~800 samples, 1645 variants) in seven scenarios
and checks that each one produces well-formed summary statistics.

It is a *smoke* test. It answers "does every code path still run and emit the
expected schema", not "are the estimates right".

```bash
tests/smoke_test.sh --list                       # show the matrix
tests/smoke_test.sh                              # all 7 scenarios (~15 min)
tests/smoke_test.sh --only vcf_glm --keep
tests/smoke_test.sh --container ghcr.io/gp2code/progwas-pipeline:1.0.0
```

Needs `nextflow`, Docker, `plink2` on `PATH` (to build the `.bed` fixture), and a
populated `REFERENCE_DIR` (hg19 + hg38 FASTA and the hg19→hg38 chain).

Each scenario writes its own params YAML under `<outdir>/params/<id>.yml`, so a
failing scenario can be reproduced by hand with a single command:

```bash
nextflow run main.nf -profile standard -params-file <outdir>/params/<id>.yml
```

## The matrix

| Scenario | Input | Model | Interaction | Genotype IDs | Pheno/covar IDs |
| --- | --- | --- | --- | --- | --- |
| `vcf_glm` | VCF | GLM linear + logistic | off | psam IID | covar `#FID`, pheno IID |
| `vcf_cph_interact` | VCF | CoxPH | on | psam IID | covar IID, pheno `#FID` |
| `vcf_gallop` | VCF | GALLOP | set, must be ignored | psam IID | covar `#FID`, pheno IID |
| `plink_glm_interact` | pgen | GLM linear + logistic | on | psam IID | covar `#FID`, pheno IID |
| `plink_fid_gallop` | pgen | GALLOP | n/a | psam FID + IID | covar `#FID`, pheno IID |
| `bed_cph` | bed | CoxPH | off | fam FID + IID | covar `#FID`, pheno IID |
| `skip_pop_split_glm` | VCF | GLM linear + logistic | off | psam IID | covar `#FID`, pheno IID |

`skip_pop_split_glm` is the only run with `skip_pop_split: true`
(`LD_PRUNE_CHR` + `SIMPLE_QC` instead of `GWASQC`). It is the branch the GP2
configs use in production, so it is in the core set even though the example
genotypes are not really ancestry-specific — that flag just skips ancestry
inference.

## Why seven runs and not forty-eight

The naive matrix is `3 inputs × 4 models × 2 interaction × 2 ID shapes = 48`
(32 if you count only VCF and PLINK and treat `.bed` as part of the PLINK path).
Most of those cells are not independent:

| Axis | Full size | Effective size | Why |
| --- | --- | --- | --- |
| Input format | 3 (vcf / pgen / bed) | pairwise with model | Format only changes the genetic-QC front end (`SPLIT_VCF`/`GENETICQC`/`MERGER_CHUNKS` vs `GENETICQCPLINK`). Everything downstream of the `chrsqced` channel consumes identical pgen, so format × model is not a real interaction. `bed` is a sub-variant of the PLINK path — it differs only in `--bfile` vs `--pfile` inside `bin/process_plink.sh`. |
| Model | 4 (linear / logistic / gallop / cph) | 3 | Logistic is **not** a separate model. `plink2 --glm` picks linear vs logistic per phenotype, so one GLM run with `pheno_name: "y,y2"` (y continuous, y2 binary 0/1) exercises both at zero extra cost. |
| Interaction | 2 | 2, GLM/CoxPH only | `modules/gwas.nf` never passes `covar_interact` to `gallop.py`, so "GALLOP × interaction" is not a cell. It is worth one *negative* assertion that the parameter is ignored rather than half-applied. |
| Phenotype/covariate FID | 2 | **0 — no longer an axis** | See below. |
| Genotype FID | 2 (3 counting `.fam`) | 3, one run each | Resolved in the genetic-QC front end by `bin/normalize_psam_iid_only.sh`. A `.psam` may be IID-only or carry FID, and a PLINK1 `.fam` always carries FID. |

### Phenotype/covariate FID is no longer a test axis

Every reader of a user-supplied phenotype or covariate file now normalises
sample identifiers through one shared implementation:

- `bin/sample_ids.py` → `normalize_sample_ids()`, used by
  `make_analysis_sets.py`, `export_plink_preprocess.py`, `make_tableone.py` and
  `gallop.py`
- `fix_plink_headers()` in `bin/survival.R`, the R counterpart

Both accept `IID` or `#IID` as the identifier, drop any `FID`/`#FID` column, and
fail loudly when no identifier column is present. `#FID IID …` and `IID …` are
therefore behaviourally identical inputs, so the shape no longer needs covering
once per model — it rides along on scenarios chosen for other reasons.

The suite still covers all four combinations at zero extra cost: the shipped
`example/covariates.tsv` has `#FID IID`, the shipped phenotype files are
IID-only, and `vcf_cph_interact` swaps both (IID-only covariates, FID-bearing
phenotypes, with `FID != IID` so a failure to drop it would be visible).

This is only true for phenotype and covariate files. Genotype FID still matters
and is handled separately — `normalize_psam_iid_only.sh` strips it during
genetic QC, and `bin/make_keep_iid.py` reconstructs it from the `.psam` when
plink2 needs a two-column `--keep` file.

### Runtime

`genetic_data_id` is pinned per *input shape* rather than per scenario, so
scenarios sharing an input reuse the cached genetic QC under
`<project_dir>/genotypes/<id>/chromosomes/`. The three VCF scenarios share one
cache; only the first pays for genetic QC.

Observed per-scenario wall time on a laptop (Docker, 2 cpus per task): GLM and
CoxPH scenarios 75–110 s each, GALLOP scenarios ~155 s each.

## What is asserted

Per scenario:

- `nextflow` exits 0.
- Every requested phenotype produced at least one `*_allresults.tsv`.
- Each result file has rows, and carries the harmonised schema shared by all
  three models: `#CHROM POS ID A1 BETA SE P P_INT P_2DF INTERACTION MODEL`.
- No duplicated column names in the header.
- At least one non-NA `P`.
- At least half of the variants produced a `P` value — a file full of
  `ERRCODE=VIF_TOO_HIGH` rows otherwise looks perfectly well-formed.
- Interaction requested → `P_INT` populated and `INTERACTION` equals the
  covariate name. Not requested → `P_INT` stays NA.
- GALLOP → `P_INT` populated and `INTERACTION` equals `TIME` (see below).
- GLM with `y,y2` → the split results carry `T_STAT`/`T_OR_F_STAT` for `y` and
  `Z_STAT`/`Z_OR_F_STAT` for `y2`, i.e. one run really did fit both a linear and
  a logistic model.
- `prepared_data` contains a `*_filtered.pca.harmonized.tsv`.

Note on GALLOP: the `*_INT` columns do **not** mean `covar_interact` there.
GALLOP fits SNP effects on both the intercept and the slope, and reuses
`BETA_INT / SE_INT / P_INT / CORR_INT / P_2DF` for the SNP×TIME slope term, with
`INTERACTION` set to the literal `TIME`. That is why `vcf_gallop` sets
`covar_interact` and then asserts `INTERACTION == TIME`: it is the check that the
parameter was ignored outright rather than leaking into the slope term.

## Bugs this suite found (now fixed)

Both were in the `GWASGLM` interaction branch of `modules/gwas.nf`.

**1. The wrong interaction term was tested, silently.** The old code did:

```bash
INTERACTION_IDX=$((N_COVAR + 2))
plink2 ... --parameters 1-${INTERACTION_IDX} --tests 1,${INTERACTION_IDX}
INTERACT_TEST="ADDx${params.covar_interact}"
```

plink2 numbers `--glm` predictors as `1` = ADD, `2 … N+1` = the covariates,
`N+2 … 2N+1` = the `ADDx<covariate>` interaction terms — and it orders the
covariates by their **column order in the covariate file**, not by the order
given to `--covar-name`. So `N_COVAR + 2` is always the interaction with
whichever covariate happens to come first, while the awk post-processing looked
for the row labelled `ADDx<covar_interact>`.

Unless `covar_interact` was the first covariate column, the run silently produced
`BETA_INT / SE_INT / P_INT` all `NA` plus a populated `P_2DF` that jointly tested
`ADD` and the **wrong** interaction term.

Verified against plink2 directly (covariates `C1..C4` in file order, target
`C2`): `--covar-name C4,C3,C2,C1` still emitted `ADD, C1, C2, C3, C4, ADDxC1 …`,
confirming file-order numbering.

The fix computes the interaction covariate's rank in file order and selects only
that term:

```bash
INTERACTION_IDX=$((N_COVAR + 1 + INTERACT_RANK))
--parameters 1-$((N_COVAR + 1)),${INTERACTION_IDX}
--tests 1,$((N_COVAR + 2))
```

`--tests` indices refer to positions *after* `--parameters` filtering — verified
empirically: `--parameters 1-5,7 --tests 1,7` yields `USER_1DF` (wrong) while
`--tests 1,6` yields `USER_2DF` (correct). A `covar_interact` that is not among
the model covariates now aborts the task instead of yielding NA columns.

`conf/examples/test_cs_linear_interact.yml` passed only because it used
`covar_interact: "SEX"`, which happened to be first. The smoke test deliberately
uses `age_at_baseline` (second in file order) so this stays covered.

**2. Duplicated output columns.** The interaction awk wrote
`<out>.<pheno>.glm.linear`, and the subsequent rename loop then matched that same
file and appended the seven `BETA_INT … MODEL` columns a second time. The padding
step is now skipped when the interaction branch already added them.

**3. The interaction covariate was left uncentered, costing ~98% of variants.**
Fixing (1) exposed this: with the correct term finally being tested, 1525 of 1558
variants came back `ERRCODE=VIF_TOO_HIGH` and no estimate at all.

`bin/export_plink_preprocess.py` deliberately skips standardising the interaction
covariate, so the interaction coefficient stays in the covariate's original units.
But it also left it *uncentered*. `plink2 --glm interaction` builds an explicit
SNP × covar column, and when the covariate has a large non-zero mean
(`age_at_baseline`, mean 72.9) that column is nearly a multiple of the SNP column.
The resulting collinearity trips plink2's variance-inflation gate (default
`--vif 50`) and the variant is dropped rather than analysed.

The fix centres the interaction covariate without scaling it. The interaction
estimate is unchanged and keeps its original units; the SNP main effect is simply
redefined as the effect at the mean covariate value, which is the usual
parameterisation for an interaction model anyway. After the fix all 1558 variants
are analysed and no row carries an ERRCODE.

Note that (1) and (3) masked each other: the old code always tested the
interaction with whichever covariate came first, which in the shipped example was
`SEX` — binary, so centering barely mattered and the VIF gate never fired.

CoxPH interaction (`bin/survival.R`) builds the model formula explicitly, has no
VIF gate, and was affected by none of the three.

### Known cosmetic difference

In interaction mode plink2 names the statistic column `T_OR_F_STAT` /
`Z_OR_F_STAT` instead of `T_STAT` / `Z_STAT`. The harmonised `BETA`/`SE`/`P`
columns are unaffected, but anything keying on the statistic column name will see
a different header between interaction and non-interaction runs. The suite accepts
either spelling.
