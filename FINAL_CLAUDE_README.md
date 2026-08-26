# pro-gwas-pipeline — Claude session handoff report

Read this before touching the codebase. It captures the pipeline's
architecture, every substantive change made across this session, what's
still open, and non-obvious Nextflow/plink2 behavior that's easy to get
wrong without re-deriving it from scratch. Companion file: `Claude_report.md`
(gitignored, local-only) has more verbose blow-by-blow detail on items 1-25
below if you need it — this file is the consolidated, portable version.

## Repo/environment context

- Local path: `/home/oiher/longGWAS/pro-gwas-pipeline`. GitHub:
  `https://github.com/Oiher/pro-gwas-pipeline`.
- **Nothing auto-commits.** The user commits/pushes manually after
  reviewing changes. Don't run `git add`/`commit`/`push` unless explicitly
  asked. Do use `git log`/`git status` (read-only) freely.
- The user runs the pipeline for real on Verily Workbench (VWB): GCS bucket
  `gs://ws-files-wb-caring-pineapple-1241`, Seqera Platform org
  `ucl-collaboration`, workspace `long-gwas`. Typical invocation loops over
  cohorts with `wb nextflow run main.nf -profile gcb_final ...`,
  `STORE_ROOT`/`PROJECT_NAME` sourced from `~/.env`, and
  `--input`/`--covarfile`/`--phenofile`/`--analysis_name` overridden per
  iteration on the command line (not baked into the yml — `-params-file`
  YAML doesn't expand env vars, see Gotchas).
- Nextflow + Docker are available in this dev sandbox. The pipeline's own
  container (`ghcr.io/hirotaka-i/long-gwas-pipeline:latest`) is pulled
  locally and has been used throughout this session to verify real
  behavior (`plink2`, `metal`, the Python/bash scripts in `bin/`) against
  synthetic fixtures, rather than trusting documentation or memory. Keep
  doing this for anything involving `plink2`/`metal` semantics — several
  real bugs this session were only caught this way (see Gotchas).

## Pipeline architecture

Nextflow DSL2 GWAS pipeline, three analysis modes: cross-sectional (GLM),
longitudinal (GALLOP/LMM), survival (Cox PH). Input is VCF or PLINK2
(pgen/pvar/psam) genotypes, plus covariate + phenotype TSVs.

```
CHECK_REFERENCES
  → genetic QC per chromosome/chunk:
      VCF path:   SPLIT_VCF → GENETICQC → MERGER_CHUNKS
      PLINK path: GENETICQCPLINK
  → MERGER_CHRS (genome-wide merge)
  → sample-level QC (one of two paths, see below):
      GWASQC (ancestry inference)  or  SIMPLE_QC (skip_pop_split: true)
  → MAKEANALYSISSETS / COMPUTE_PCA / HARMONIZE_CATEGORICAL_COVARS
      (analysis-specific prep — first point covarfile is consulted)
  → EXPORT_PLINK (cross-sectional) or RAWFILE_EXPORT (longitudinal/survival)
      (first point phenofile is consulted — samples without a phenofile
      row are dropped here, via an inner join, not earlier)
  → GWASGLM / GWASGALLOP / GWASCPH (the actual regression)
  → SAVEGWAS / MANHATTAN / TABLEONE
```

**Important**: variant-level QC, merging, and sample-level QC (ancestry/
kinship/PCA — the most expensive stage) all run on the **full genotype
input sample count**, not the subset that actually has phenotype data.
Phenotype-file filtering doesn't happen until `EXPORT_PLINK`/
`RAWFILE_EXPORT`, the very last prep step. If `--input` has 11k samples but
`--phenofile` only covers 2k, compute cost for everything before that last
step scales with 11k, not 2k. This is architecturally sound (kinship/
ancestry/PCA estimates benefit from a larger reference pool), just worth
knowing when estimating run cost.

`run_metal.nf` and `run_focus.nf` are separate, standalone entrypoints
(their own `nextflow.enable.dsl=2`, not wired into `main.nf`'s DAG) —
`run_metal.nf` for cross-cohort/cross-population meta-analysis of already-
produced `*_allresults.tsv` files (wraps the `metal` CLI), `run_focus.nf`
for targeted re-analysis of a small variant set without re-running full
QC/PCA.

## Caching model (correctness-relevant, not just performance)

- `genotypes/${genetic_cache_key}/chromosomes/` and
  `analyses/${genetic_cache_key}/genetic_qc/` are shared across every run
  using the same `genetic_cache_key` (auto-computed from
  format/ancestry/assembly/maf/kinship, or overridden via
  `--genetic_data_id`).
- `analyses/${genetic_cache_key}/${analysis_name}/` (prepared_data,
  gwas_results) is always `analysis_name`-scoped — recomputed per analysis
  regardless of cache-key sharing.
- **`genetic_data_id` must be unique per distinct *set of input files*, not
  per abstract "cohort"** — reusing an ID across two runs with genuinely
  different `--input` file sets previously caused silent cross-
  contamination (fixed, item 6 below; a `log.warn` now catches this).
- Reusing `genetic_data_id` across several *different-phenotype* analyses
  on the *same* genotype input (e.g. one GLM run on UPDRS, one CoxPH run on
  SEADL, one LMM run on MoCA, same 11k-sample input) is the intended way to
  share expensive genetic QC. But **not every shared-cache-key stage skips
  automatically**:
  - `GENETICQCPLINK`/`GENETICQC` use `storeDir` → unconditional skip-if-
    exists, reused regardless of `-resume`.
  - `MERGER_CHRS`, `GWASQC`, `SIMPLE_QC` use `publishDir` + `cache 'deep'`
    only → reused **only if `-resume` is passed** on the later runs, since
    this relies on Nextflow's normal task-hash cache, not an unconditional
    file check. Without `-resume`, these re-execute (correct results,
    wasted compute).
  - `MAKEANALYSISSETS` onward is `analysis_name`-scoped and genuinely
    re-executes per analysis (correct — different phenofile/covarfile
    content).
  - **Practical recommendation**: same `genetic_data_id`, distinct
    `--analysis_name` per analysis, and pass `-resume` on every run after
    the first to actually get the `MERGER_CHRS`/`GWASQC` savings.

## QC steps applied to input genotypes

Four stages, each on a different unit of data (kinship/PCA/ancestry
calculation-support filters — e.g. the LD-pruning used only to select SNPs
for computing a metric — are called out separately, not counted as
"dataset QC filters" below):

**1. Variant-level** (`bin/process1.sh` for VCF input, `bin/process_plink.sh`
for PLINK input — parallel scripts, note their positional CLI args are in
a *different order* between the two, a real pre-existing landmine if you
ever add a new positional arg to either): PASS/imputation-R2 filter,
liftover to hg38, normalize/split-multiallelic, biallelic-SNP-only +
`--mac 2` + `--geno 0.1` + REF/ALT-aligned + dedup.

**2. Merge** (`MERGER_CHUNKS` → `MERGER_CHRS`): no filtering, structural
only.

**3. Sample-level**, one of two paths depending on `skip_pop_split`:
- `GWASQC` (`bin/addi_qc_pipeline.py`/`bin/qc.py`, ancestry-inference
  path): `--mind 0.05` call-rate, heterozygosity outliers at absolute
  `F ≤ -0.25` or `F ≥ 0.25`, ancestry projection against a reference panel,
  kinship via `--king-cutoff 0.0884` (**hardcoded — `params.kinship` does
  NOT apply on this path**, only on `SIMPLE_QC` below; worth fixing for
  consistency if it matters to you).
- `SIMPLE_QC` (`bin/simple_qc.sh`, `skip_pop_split: true`, no ancestry
  inference): `--mind 0.05`, heterozygosity outliers at `±3 SD` from the
  mean (**different methodology than `GWASQC`'s absolute bound** — not
  just a hardcoded-value difference, a genuinely different formula),
  kinship via `--king-cutoff ${params.kinship}` (this is the only path
  where the param actually takes effect).

**4. Final export filters** (`EXPORT_PLINK`/`RAWFILE_EXPORT` prep, actual
filtering happens in the `plink2 --glm` calls in `modules/gwas.nf` for
cross-sectional, or in `RAWFILE_EXPORT` itself for longitudinal/survival):
`--maf ${params.minor_allele_freq}`, `--mac ${params.minor_allele_ct}`,
`--hwe 1e-6` (hardcoded).

### Customizable vs. hardcoded (excluding kinship/PCA/ancestry machinery)

| Customizable | Default | Hardcoded | Value |
|---|---|---|---|
| `r2thres` | `-9` (disabled) | PASS filter | `.,PASS` |
| `minor_allele_freq` | `0.05` | `--mac 2` (singleton) | variant-level |
| `minor_allele_ct` | `20` | `--geno 0.1` | variant missingness |
| | | `--mind 0.05` | sample call rate (both paths) |
| | | heterozygosity bound | `±0.25` abs. (`GWASQC`) / `±3 SD` (`SIMPLE_QC`) |
| | | `--hwe 1e-6` | final export |

**`-9` does NOT mean "disabled" for `--maf`/`--mac`/`--hwe`** — verified
directly against the real `plink2` binary: `--maf -9`/`--mac -9` hard-error
(`must be >= 0`) immediately. Only `r2thres` has genuine `-9`-means-disabled
behavior, and only because `process1.sh`/`process_plink.sh` have explicit
bash conditionals checking `r2 > 0` before applying the filter at all — it's
not a `plink2` convention. **The correct no-op value for `--maf`/`--mac`/
`--hwe` is `0`** — also verified directly: these are lower-bound threshold
filters (exclude values *below* the threshold), so `0` genuinely keeps
everything (confirmed with synthetic fixtures including a variant with a
severe HWE violation, p=1.34e-6 — `--hwe 0.001` correctly removed it,
`--hwe 0` correctly kept it).

There's an unimplemented **plan** to expose `hwe`/`mind`/`geno`/
heterozygosity-bound as real pipeline params — see `/home/oiher/.claude/plans/i-will-take-note-typed-hejlsberg.md`
for the full per-parameter breakdown (exact call sites, a real
`storeDir`-staleness risk specific to `geno`'s `GENETICQCPLINK` call site,
and an open design decision on whether to unify `GWASQC`/`SIMPLE_QC`'s two
different heterozygosity formulas or parameterize them independently).
Not started — planning only, confirmed with the user.

## Changes made this session (summary — see `Claude_report.md` for full detail on 1-25)

**GCB/Google Batch resourcing (items 1-16):**
1. Unified `gcb.config`+`gcb2.config` → `conf/profiles/gcb_final.config`.
2. Kept `singularity{}` block (inert under `google-batch`, documented why).
3. Full disk-sizing overhaul across `small`/`medium`/`large`/`very_large`
   labels + per-process `withName` overrides, driven by a real
   `SSD_TOTAL_GB` quota hit (500GB, `europe-west4`).
4. Upfront `--input`/`--covarfile`/`--phenofile` null-param validation.
5. Empty-input-glob guard (`.ifEmpty{error(...)}`).
6. **Cache-contamination guard** — the biggest correctness fix: prevents
   silently merging unrelated samples when a `genetic_data_id` is reused
   across different `--input` file sets.
7. Replaced fragile bare-env-var Nextflow pattern with explicit
   `params.analyses_dir`/`params.genotypes_dir`.
8. `GENETICQCPLINK` right-sizing — also fixed a `two_cpu_large_mem` label
   with a physically-unsatisfiable `cpus=2,memory=128GB` ratio (real GCP
   ceiling is ~8GB/vCPU).
9. Phenotype/covariate column-presence validation in `main.nf` (the
   ancestor of items 17-19 below).
10. `maxForks` gap audit across every per-chromosome/per-chunk fan-out
    process.
11. Manhattan plot variant labeling — added, iterated through several
    `qmplot` library gotchas, ultimately **removed entirely** after label-
    to-point correspondence proved unreliable on real production data. Plot
    now only draws `NOMINAL_P`(1e-5)/`MTC_P`(5e-8) reference lines.
12. `GWASGLM` `manifest.tsv` publish collision (22 chromosome tasks writing
    the same filename) — fixed with per-chromosome unique naming; also hit
    and fixed a `def`-scoping-invisible-to-`output:` Nextflow gotcha along
    the way (see Gotchas).
13. Cross-process combined disk budget — `RAWFILE_EXPORT` pipelines
    directly into `GWASGALLOP`/`GWASCPH`, so their `maxForks` caps must be
    budgeted together, not independently.
14. User-facing `CODE_GCE_QUOTA_EXCEEDED` explanation in `main.nf`/README.
15. `gcb_scaleable` profile + `bin/detect_gcp_quota.sh` — makes the
    500GB/200CPU assumption configurable/auto-detectable for other GCP
    projects.
16. Automatic quota detection printed at every `main.nf` launch (superseded
    item 15's static note).

**Reactive fixes from real multi-phenotype production runs (items 17-25):**
17. Survival phenotype validation — require `tstart`/`tend` columns +
    binary `pheno_name` coding when `survival_flag`.
18. `time_col` requirement scoping fix (mandatory for `longitudinal_flag`
    only, not `survival_flag`).
19. **`--pheno_name` comma-vs-space delimiter consolidation** — one shared
    `parsePhenoNames()` helper, fixed 4 independently-drifted call sites.
20. `TABLEONE` per-phenotype fan-out (was silently incomplete for
    multi-phenotype runs).
21. `GWASGLM` input staging guard (`awk` field-count self-check before
    `plink2` runs).
22. `export_plink_preprocess.py` self-check for malformed output.
23. **Real root cause of a recurring `plink2` "fewer tokens" error**:
    missing quantitative phenotype/covariate values were left as empty
    strings (pandas default), which `plink2`'s tokenizer collapses instead
    of treating as valid fields — fixed by filling with `plink2`'s `NA`
    token instead.
24. `run_metal.nf`: cross-cohort meta-analysis via new `--metal_pheno_name`
    param (grouped-by-phenotype fan-out, one METAL run per phenotype).
25. `run_metal.nf`: identical-basename collision fix — every cohort's
    output file shares the same basename (cohort identity lives only in
    the parent directory), fixed via `stageAs: 'input??/*'` numbered
    subdirectory staging.

**26. `export_plink_preprocess.py`: all-missing-phenotype miscoding bug**
(this session, not yet in `Claude_report.md`'s numbered list). A
phenotype with **zero non-missing values** for a given study-arm export
(e.g. `UPDRS_pI` entirely `NaN` for a specific small site) hit a Python
logic gap: `set(unique_vals).issubset({1, 2, -9})` is vacuously `True` for
an empty set, and that branch (unlike the sibling binary-0/1 branch just
above it) had no `len(unique_vals) > 0` guard — so an all-missing
phenotype fell into the "already PLINK format" branch and got every value
coded as literal `-9` instead of falling through to the (already-existing)
"fill with `NA`" branch. Fixed by adding the missing `len() > 0` guard.
Verified directly against real `plink2` — reproduced the user's exact
error (`All samples for --glm phenotype 'X' are controls`) with a minimal
fixture, confirmed the fix produces `NA` instead of `-9`.

**Important, still-open finding from the same investigation**: fixing the
encoding bug alone does **not** fully resolve this class of failure.
Verified separately: `plink2`'s multi-phenotype `--glm` call aborts
**entirely** if *any* requested phenotype has zero valid values — this
happens regardless of whether the missing values are coded as `NA` or
`-9`, so it's not an encoding issue at all. Dropping the empty phenotype
from `--pheno-name` for that specific study-arm/chromosome task is what
actually lets `plink2` proceed. The pipeline currently builds one static
`--pheno-name` list (from `params.pheno_name`) applied identically to
every chromosome × study-arm task — it has no mechanism to detect "this
phenotype has zero data for this particular study arm" and drop it from
that task's list. See Pending below.

## Pending / could be addressed in future sessions

- **`plink2` multi-phenotype all-missing abort** (new, see item 26 above,
  not yet resolved or planned in detail). Would need `export_
  plink_preprocess.py` (which already computes per-phenotype valid-value
  counts) to emit a per-study-arm "phenotypes with data" list, and
  `GWASGLM` (`modules/gwas.nf`) to build its `--pheno-name` from that
  instead of the static pipeline-wide list. Discussed with the user but
  not planned/scoped — they may instead choose to curate which
  cohorts/sites go into which phenotype's analysis manually.
- **QC threshold customization plan** (see the QC section above and the
  saved plan at `/home/oiher/.claude/plans/i-will-take-note-typed-hejlsberg.md`)
  — planned, not implemented. Exposes `hwe`, `mind`, `geno`, and
  heterozygosity bounds as real params. Has one real open design decision
  (unify vs. independently parameterize the two different heterozygosity
  formulas) and one real risk to handle carefully (`geno`'s
  `GENETICQCPLINK` call site uses `storeDir`, so changing it silently
  reuses stale cached chromosomes unless `genetic_data_id` also changes).
- **`network`/`subnetwork` hardcoding in `gcb_final.config`** — unclear if
  a guaranteed VWB convention or project-specific; explicitly deferred,
  needs real GCP access to check across workspaces.
- **`GENETICQC`'s unused `fOrig` staging** — full per-chromosome VCF staged
  into every chunk task but never read; dropping it would cut real disk
  need from ~60GB toward ~15-20GB on the highest-concurrency process in
  the pipeline. Deferred as a follow-up, not implemented.
- **`two_cpu_large_mem` label's unsatisfiable `128GB/2cpu` ratio** — only
  fixed for `GENETICQCPLINK`. `MAKEANALYSISSETS`/`SAVEGWAS`/`MANHATTAN`
  still carry it, masked so far by low concurrency.
- **GCP quota increase** — the actual permanent fix for `SSD_TOTAL_GB`;
  not something Claude can do, needs the user/a project admin.
- **`adwb.config` label mismatch** (`large_mem` vs. `large`/
  `two_cpu_large_mem`) — pre-existing, unrelated to this session, not
  fixed.
- **Multi-phenotype consolidation recommendation** — user currently runs
  one `--input`/`genetic_data_id` per phenotype due to non-overlapping
  per-phenotype genotype exports; recommended consolidating to one
  genotype file + one multi-column phenofile per cohort instead, to share
  genetic QC across phenotypes. Not adopted — requires the user to change
  their upstream data-prep, not just this repo.
- **Manhattan/QQ plots for `run_metal.nf`'s meta-analysis output** — `.TBL`
  files use different column names/naming than `bin/manhattan.py` expects
  (`MarkerName` vs `#CHROM`/`POS`, `P-value` vs `P`). User explicitly asked
  to remember this for a future session, not implement now.
- **`GWASQC`'s hardcoded `--king-cutoff 0.0884`** ignores `params.kinship`
  entirely (only `SIMPLE_QC` respects it) — noted during the QC-parameter
  audit, not yet decided whether this inconsistency should be fixed.

## Nextflow/plink2 gotchas (save time on future changes)

- **`Channel.combine()` auto-flattens a bare `List`** into separate
  positional tuple elements. Wrap it one level deeper (`.map{ tags ->
  [tags] }`) before combining if you need it to stay one element.
- **`Channel.join()` collapses multiple left-side entries sharing a key**
  down to one match, rather than pairing every left entry with the right
  match. Use `combine()` + `filter{}` instead if you need to preserve all
  of the left side's entries.
- **`withName` overrides `withLabel` per-directive**, not wholesale —
  unset directives in a `withName` block still fall through to the
  matching `withLabel` block.
- **`-params-file` YAML does not resolve environment variables** — only a
  small set of Nextflow's own built-ins (`${projectDir}` etc.). Pass
  env-var-derived paths via `--param` on the CLI instead, where the shell
  expands them first.
- **`nextflow.config`'s DSL sandbox rejects user-defined functions doing
  arithmetic** (`Unknown method invocation 'multiply' on Integer type`)
  even though identical arithmetic works as a top-level statement. Use
  inline top-level `def` computations instead of a shared helper function
  inside `.config` files specifically (this restriction does *not* apply
  inside `.nf` scripts, which run with the full Groovy runtime — confirmed
  by testing both contexts).
- **A `def`-scoped local declared in a process `script:` block is invisible
  to that same process's `output:` block** — `output:` silently evaluates
  it as the string `"null"` rather than erroring. Don't use `def` for any
  variable an `output:` path pattern references.
- **`waitForOrKill()`'s return value is `null` unconditionally** — on both
  a killed-by-timeout process and one that completed normally. Check
  `exitValue()` instead (always safely readable after `waitForOrKill()`
  returns, even when killed).
- **GCS paths are case-sensitive**; the VWB-mounted view
  (`/home/jupyter/workspace/ws_files/...`) mirrors the bucket 1:1 including
  case, but Google Batch tasks need the real `gs://` URI, not the mount
  path (the mount only exists on the interactive VM).
- **GCP machine types cap memory-per-vCPU at ~8GB** (even n2-highmem) — a
  directive like `cpus=2, memory='128 GB'` silently provisions a much
  bigger machine regardless of the stated `cpus`, actively misleading both
  cost and `CPUS`-quota accounting. Check the ratio against real GCP
  families (~4GB/vCPU standard, ~8GB/vCPU highmem, ~1-2GB/vCPU highcpu)
  whenever sizing a `gcb_final` label.
- **`plink2 --maf`/`--mac`/`--hwe` are lower-bound filters, not exact-match
  or exclusion filters** — `0` is the genuine no-op (keeps everything);
  `-9` is not a recognized "disabled" sentinel for these three and hard-
  errors immediately (`must be >= 0`). Only `r2thres` has real `-9`-means-
  disabled behavior, implemented as pipeline-side bash conditionals, not a
  `plink2` convention.
- **`plink2`'s multi-phenotype `--glm` aborts entirely if any one requested
  phenotype has zero valid values** — not a per-phenotype skip, the whole
  command fails, taking down every other (valid) phenotype in the same
  call too. See item 26/Pending above.
- **An empty Python set is vacuously a subset of any set**
  (`set().issubset(X)` is always `True`) — a real bug source when checking
  "are all values in the allowed set" without also checking there *are*
  values (`len(unique_vals) > 0`). Bit `export_plink_preprocess.py` once
  already (item 26); worth grepping for the same pattern elsewhere if
  similar all-missing-column bugs turn up.
- **Identical basenames across different source directories** (e.g. every
  cohort's `EUR_case_UPDRS_pI_allresults.tsv`) will collide when Nextflow
  tries to stage them flat into one task directory. Use `path(files,
  stageAs: 'input??/*')` to number them into separate subdirectories.
