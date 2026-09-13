# Pipeline smoke test

`tests/smoke_test.sh` runs the main pipeline (`main.nf`) end to end over the
bundled example data (chr20–22, ~800 samples, 1645 variants) in a handful of
scenarios and checks that each one produces well-formed summary statistics.

It is a *smoke* test. It answers "does every code path still run and emit the
expected schema", not "are the estimates right".

```bash
tests/smoke_test.sh --list                       # show the matrix
tests/smoke_test.sh                              # 6 core scenarios (~10 min)
tests/smoke_test.sh --tier all                   # + bed and skip_pop_split
tests/smoke_test.sh --only vcf_glm --keep
tests/smoke_test.sh --container ghcr.io/gp2code/progwas-pipeline:1.0.0
```

Needs `nextflow`, Docker, and a populated `REFERENCE_DIR` (hg19 + hg38 FASTA and
the hg19→hg38 chain). `plink2` on `PATH` is only needed for the optional
`bed_glm` scenario.

Each scenario writes its own params YAML under `<outdir>/params/<id>.yml`, so a
failing scenario can be reproduced by hand with a single command:

```bash
nextflow run main.nf -profile standard -params-file <outdir>/params/<id>.yml
```

## Why the matrix is only six runs

The naive matrix is `2 inputs × 4 models × 2 interaction × 2 ID shapes = 32`.
Most of those cells are not independent:

| Axis | Full size | Effective size | Why |
| --- | --- | --- | --- |
| Input format | 2–3 (vcf / pgen / bed) | 2 core, 3 with `--tier all` | Format only changes the genetic-QC front end (`SPLIT_VCF`/`GENETICQC`/`MERGER_CHUNKS` vs `GENETICQCPLINK`). Everything downstream of the `chrsqced` channel consumes identical pgen, so format × model is not a real interaction — pairwise coverage suffices. |
| Model | 4 (linear / logistic / gallop / cph) | 3 | Logistic is **not** a separate model. `plink2 --glm` picks linear vs logistic per phenotype, so one GLM run with `pheno_name: "y,y2"` (y continuous, y2 binary 0/1) exercises both at zero extra cost. |
| Interaction | 2 | 2, but GLM/CoxPH only | `modules/gwas.nf` never passes `covar_interact` to `gallop.py`, so "GALLOP × interaction" is not a cell. It is worth one *negative* assertion that the parameter is ignored rather than half-applied. |
| Sample IDs | 2 | 2, covered once each | Both genetic-QC paths call `bin/normalize_psam_iid_only.sh`, so FID is stripped before anything model-specific sees it. The axis needs covering once on the genotype side (an FID-bearing `.psam`) and once on the covariate side (a covariate file with no `#FID`) — not once per model. |

The six core scenarios form a pairwise covering array over format × model,
model × interaction, and both ID shapes:

| Scenario | Input | Model | Interaction | psam IDs | covar IDs |
| --- | --- | --- | --- | --- | --- |
| `vcf_glm` | VCF | GLM linear + logistic | off | IID | `#FID` |
| `vcf_cph_interact` | VCF | CoxPH | on | IID | IID only |
| `vcf_gallop` | VCF | GALLOP | set, must be ignored | IID | `#FID` |
| `plink_fid_glm_interact` | pgen | GLM linear + logistic | on | FID + IID | IID only |
| `plink_cph` | pgen | CoxPH | off | IID | `#FID` |
| `plink_gallop` | pgen | GALLOP | n/a | IID | IID only |

Optional (`--tier all`):

| Scenario | Covers |
| --- | --- |
| `bed_glm` | the `.bed/.bim/.fam` branch of `bin/process_plink.sh` |
| `skip_pop_split_glm` | `skip_pop_split: true` — `LD_PRUNE_CHR` + `SIMPLE_QC` instead of `GWASQC` |

Both extras were verified to pass on the example data. `skip_pop_split: true` in
particular does run here, even though the example genotypes are not really
ancestry-specific — it simply skips ancestry inference. That branch is a large
one (`LD_PRUNE_CHR`, `MERGER_CHRS`, `SIMPLE_QC` instead of `GWASQC`) and is what
the GP2 configs actually use in production, so it is worth running periodically
even though it is not in the core tier.

The `bed_glm` fixture is built with `plink2 --max-alleles 2`, because the example
pgen contains multiallelic variants that a PLINK1 `.bim` cannot represent.

Runtime is kept down by pinning `genetic_data_id` per *input shape* rather than
per scenario, so scenarios that share an input reuse the cached genetic QC under
`<project_dir>/genotypes/<id>/chromosomes/`. The first VCF scenario pays for
genetic QC; the other two VCF scenarios skip it.

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
- Interaction requested → `P_INT` populated and `INTERACTION` equals the
  covariate name. Not requested → `P_INT` stays NA.
- GALLOP → `P_INT` populated and `INTERACTION` equals `TIME` (see below).
- GLM with `y,y2` → the split results carry `T_STAT`/`T_OR_F_STAT` for `y` and
  `Z_STAT`/`Z_OR_F_STAT` for `y2`, i.e. one run really did fit both a linear and
  a logistic model.
- `prepared_data` contains a `*_filtered.pca.harmonized.tsv`.

Note on GALLOP: the `*_INT` columns do **not** mean `covar_interact` there. GALLOP
fits SNP effects on both the intercept and the slope, and reuses `BETA_INT /
SE_INT / P_INT / CORR_INT / P_2DF` for the SNP×TIME slope term, with
`INTERACTION` set to the literal `TIME`. That is why `vcf_gallop` sets
`covar_interact` and then asserts `INTERACTION == TIME`: it is the check that the
parameter was ignored outright rather than leaking into the slope term.

## Known failure: `plink_fid_glm_interact`

This scenario fails against the current pipeline, and the failure is real, not a
test artefact. `modules/gwas.nf` (GWASGLM, interaction branch) does:

```bash
INTERACTION_IDX=$((N_COVAR + 2))
plink2 ... --parameters 1-${INTERACTION_IDX} --tests 1,${INTERACTION_IDX}
INTERACT_TEST="ADDx${params.covar_interact}"
```

`--parameters 1-(N+2)` always selects the interaction with whichever covariate
plink2 orders **first**, and plink2 orders covariates by their column order in
the covariate file, not by the order given to `--covar-name`. Meanwhile the awk
post-processing looks for the row labelled `ADDx<covar_interact>`.

So unless `covar_interact` happens to be the first covariate column in
`*_filtered.pca.pheno.tsv`, the run silently produces:

- `BETA_INT / SE_INT / P_INT` all `NA` (the `ADDx<covar_interact>` row is never
  emitted), and
- a populated `P_2DF` that is the joint test of `ADD` plus the **wrong**
  interaction term.

Reproduced directly with plink2 on the pipeline's own prepared data
(`covar_interact: age_at_baseline`, covariate file order `SEX, age_at_baseline,
PC1, …`): plink2 emitted `ADD, SEX, age_at_baseline, PC1, PC2, PC3, site_B,
site_C, specimen_no, ADDxSEX, USER_2DF` — `ADDxSEX`, not
`ADDxage_at_baseline`.

`conf/examples/test_cs_linear_interact.yml` passes only because it uses
`covar_interact: "SEX"` and `SEX` happens to be the first covariate column.

A second, smaller defect shows up in the same branch: the interaction awk writes
`<out>.<pheno>.glm.linear`, and the subsequent rename loop then matches that same
file and appends the seven `BETA_INT … MODEL` columns a second time. Interaction
results therefore carry each of those columns twice, the second copy all `NA`.

CoxPH interaction (`bin/survival.R`) builds the model formula explicitly and is
not affected — `vcf_cph_interact` passes.
