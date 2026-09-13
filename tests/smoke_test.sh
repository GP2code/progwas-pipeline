#!/usr/bin/env bash
#
# Smoke test for the main proGWAS pipeline (main.nf).
#
# Runs a small covering set of scenarios over the bundled example data
# (chr20-22, ~800 samples) and asserts that each one produces well-formed
# summary statistics. It is a *smoke* test: it checks that every code path runs
# end to end and emits the harmonised result schema. It does NOT check
# statistical correctness.
#
# Usage:
#   tests/smoke_test.sh                        # the 6 core scenarios
#   tests/smoke_test.sh --tier all             # core + optional scenarios
#   tests/smoke_test.sh --only vcf_glm,plink_gallop
#   tests/smoke_test.sh --list
#   tests/smoke_test.sh --dry-run
#
# Options:
#   --tier core|all      Which scenario set to run (default: core)
#   --only a,b,c         Run only these scenario ids (overrides --tier)
#   --profile NAME       Nextflow profile (default: standard)
#   --container IMAGE    Override process.container for every task
#   --outdir DIR         Where to put results (default: .smoke_test/<timestamp>)
#   --keep               Keep the output directory even when everything passes
#   --list               Print the scenario matrix and exit
#   --dry-run            Print the nextflow commands and exit
#
# Environment:
#   REFERENCE_DIR        Reference genomes (default: <repo>/References)
#   NXF_SYNTAX_PARSER    Set to v1 automatically if unset (needed on Nextflow 26.04+)
#
# Written for bash 3.2 (the macOS system bash), so no mapfile/associative arrays.

set -o pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT" || exit 1

# ------------------------------------------------------------------ defaults
TIER="core"
ONLY=""
PROFILE="standard"
CONTAINER=""
OUTDIR=""
KEEP=0
DRY_RUN=0
LIST_ONLY=0

while [ $# -gt 0 ]; do
  case "$1" in
    --tier)      TIER="$2"; shift 2 ;;
    --only)      ONLY="$2"; shift 2 ;;
    --profile)   PROFILE="$2"; shift 2 ;;
    --container) CONTAINER="$2"; shift 2 ;;
    --outdir)    OUTDIR="$2"; shift 2 ;;
    --keep)      KEEP=1; shift ;;
    --list)      LIST_ONLY=1; shift ;;
    --dry-run)   DRY_RUN=1; shift ;;
    -h|--help)   sed -n '3,33p' "${BASH_SOURCE[0]}" | sed 's/^#\{0,1\} \{0,1\}//'; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

: "${OUTDIR:=$REPO_ROOT/.smoke_test/$(date +%Y%m%d_%H%M%S)}"
: "${REFERENCE_DIR:=$REPO_ROOT/References}"
: "${NXF_SYNTAX_PARSER:=v1}"
export REFERENCE_DIR NXF_SYNTAX_PARSER

FIXTURES="$OUTDIR/fixtures"
PARAMSDIR="$OUTDIR/params"
LOGDIR="$OUTDIR/logs"
STORE_ROOT="$OUTDIR/store"
PROJECT_NAME="smoke"
PROJECT_DIR="$STORE_ROOT/$PROJECT_NAME"

# ------------------------------------------------------------- scenario table
#
# Why the matrix is this small:
#
#   input format {vcf, plink}
#       Only changes the genetic-QC front end (SPLIT_VCF/GENETICQC/MERGER_CHUNKS
#       vs GENETICQCPLINK). Everything downstream of the `chrsqced` channel sees
#       identical pgen, so format x model is not a real interaction -- pairwise
#       coverage is enough, not the full cross product.
#
#   model {glm, cph, gallop}
#       Logistic is NOT a separate model: plink2 --glm picks linear vs logistic
#       per phenotype, so `pheno_name: "y,y2"` (y continuous, y2 binary) covers
#       both in one run at zero extra cost.
#
#   interaction {off, on}
#       Implemented for GLM and CoxPH only. modules/gwas.nf never passes
#       covar_interact to gallop.py, so gallop x interaction is not a real cell.
#       Scenario `vcf_gallop` sets covar_interact anyway and asserts the
#       interaction columns stay empty, i.e. that it is ignored outright rather
#       than half-applied.
#
#   sample IDs {IID, FID+IID}
#       Both genetic-QC paths call normalize_psam_iid_only.sh, so FID is stripped
#       before anything model-specific sees it. The axis therefore needs covering
#       once on the genotype side (an FID-bearing .psam) and once on the
#       covariate side (a covariate file with no #FID column) -- not once per
#       model.
#
# The six core scenarios are a pairwise covering array over format x model,
# model x interaction, and both ID shapes.

SCENARIOS_CORE="vcf_glm vcf_cph_interact vcf_gallop plink_fid_glm_interact plink_cph plink_gallop"
SCENARIOS_EXTRA="bed_glm skip_pop_split_glm"

scenario_desc() {
  case "$1" in
    vcf_glm)                echo "VCF  | GLM linear+logistic | no interact | psam IID     | covar #FID     | manhattan" ;;
    vcf_cph_interact)       echo "VCF  | CoxPH               | interact    | psam IID     | covar IID-only" ;;
    vcf_gallop)             echo "VCF  | GALLOP              | ignored     | psam IID     | covar #FID" ;;
    plink_fid_glm_interact) echo "pgen | GLM linear+logistic | interact    | psam FID+IID | covar IID-only" ;;
    plink_cph)              echo "pgen | CoxPH               | no interact | psam IID     | covar #FID" ;;
    plink_gallop)           echo "pgen | GALLOP              | n/a         | psam IID     | covar IID-only" ;;
    bed_glm)                echo "bed  | GLM linear+logistic | no interact | fam FID+IID  | covar #FID     | [extra]" ;;
    skip_pop_split_glm)     echo "VCF  | GLM linear+logistic | no interact | psam IID     | covar #FID     | [extra] skip_pop_split=true" ;;
    *)                      echo "(unknown scenario)" ;;
  esac
}

if [ "$LIST_ONLY" -eq 1 ]; then
  printf '%-24s %s\n' "SCENARIO" "DESCRIPTION"
  for s in $SCENARIOS_CORE $SCENARIOS_EXTRA; do
    printf '%-24s %s\n' "$s" "$(scenario_desc "$s")"
  done
  exit 0
fi

case "$TIER" in
  core) SCENARIOS="$SCENARIOS_CORE" ;;
  all)  SCENARIOS="$SCENARIOS_CORE $SCENARIOS_EXTRA" ;;
  *)    echo "Unknown --tier '$TIER' (expected core|all)" >&2; exit 2 ;;
esac

if [ -n "$ONLY" ]; then
  SCENARIOS="$(echo "$ONLY" | tr ',' ' ')"
fi

# ---------------------------------------------------------------- preflight
fail_hard() { echo "ERROR: $*" >&2; exit 1; }

command -v nextflow >/dev/null || fail_hard "nextflow is not on PATH"

for f in Genome/hg38.fa.gz Genome/hg38.fa.gz.fai Genome/hg38.fa.gz.gzi \
         Genome/hg19.fa.gz Genome/hg19.fa.gz.fai Genome/hg19.fa.gz.gzi \
         liftOver/hg19ToHg38.over.chain.gz; do
  [ -f "$REFERENCE_DIR/$f" ] || fail_hard \
    "missing reference $REFERENCE_DIR/$f -- run: bin/download_references.sh hg19 '$REFERENCE_DIR'"
done

for f in example/genotype/chr20.vcf example/genotype_plink/chr20.pgen \
         example/covariates.tsv example/phenotype.cs2.tsv \
         example/phenotype.lt.tsv example/phenotype.surv.tsv; do
  [ -f "$REPO_ROOT/$f" ] || fail_hard "missing example fixture $f"
done

mkdir -p "$FIXTURES" "$PARAMSDIR" "$LOGDIR" "$STORE_ROOT" || exit 1

# ---------------------------------------------------------------- fixtures
#
# Derived inputs the repo does not ship: a covariate file with no #FID column,
# a PLINK set whose .psam carries a real FID, and a .bed/.bim/.fam set.

build_fixtures() {
  echo "==> building fixtures in $FIXTURES"

  # Covariate file with the #FID column removed (the IID-only input shape).
  awk -F'\t' -v OFS='\t' '
    NR==1 { for (i=1;i<=NF;i++) if ($i=="#FID" || $i=="FID") fid=i }
    { out=""; for (i=1;i<=NF;i++) { if (i==fid) continue; out=(out==""?$i:out OFS $i) } print out }
  ' example/covariates.tsv > "$FIXTURES/covariates.nofid.tsv"
  if head -1 "$FIXTURES/covariates.nofid.tsv" | grep -q 'FID'; then
    fail_hard "covariates.nofid.tsv still carries an FID column"
  fi

  # PLINK set with FID+IID in the .psam. FID is deliberately different from IID
  # so normalize_psam_iid_only.sh has something real to strip.
  mkdir -p "$FIXTURES/plink_fid"
  for chr in 20 21 22; do
    cp "example/genotype_plink/chr${chr}.pgen" "$FIXTURES/plink_fid/chr${chr}.pgen"
    cp "example/genotype_plink/chr${chr}.pvar" "$FIXTURES/plink_fid/chr${chr}.pvar"
    awk -F'\t' -v OFS='\t' '
      NR==1 { out = "#FID" OFS "IID"; for (i=2;i<=NF;i++) out = out OFS $i; print out; next }
            { out = "fam_" $1 OFS $1; for (i=2;i<=NF;i++) out = out OFS $i; print out }
    ' "example/genotype_plink/chr${chr}.psam" > "$FIXTURES/plink_fid/chr${chr}.psam"
  done

  # PLINK1 binary set, for the .bed branch of bin/process_plink.sh.
  # --max-alleles 2 is required: the example pgen carries multiallelic variants
  # and a PLINK1 .bim cannot represent those.
  case " $SCENARIOS " in
    *" bed_glm "*)
      if command -v plink2 >/dev/null; then
        mkdir -p "$FIXTURES/plink_bed"
        for chr in 20 21 22; do
          plink2 --pfile "example/genotype_plink/chr${chr}" --max-alleles 2 \
                 --make-bed --out "$FIXTURES/plink_bed/chr${chr}" \
                 > "$FIXTURES/plink_bed/chr${chr}.makebed.log" 2>&1 ||
            fail_hard "plink2 --make-bed failed for chr${chr} (see $FIXTURES/plink_bed/chr${chr}.makebed.log)"
        done
      else
        echo "    plink2 not on PATH -- skipping bed_glm" >&2
        SCENARIOS="$(echo "$SCENARIOS" | sed 's/bed_glm//')"
      fi
      ;;
  esac
}

# ------------------------------------------------------------ params writer
#
# Each scenario gets its own params YAML so a failing scenario can be re-run by
# hand with one command:
#   nextflow run main.nf -profile standard -params-file <outdir>/params/<id>.yml
#
# genetic_data_id is pinned per *input shape*, not per scenario, so scenarios
# sharing an input reuse the cached genetic QC under
# <project_dir>/genotypes/<id>/chromosomes/. That is what keeps the suite to a
# few minutes rather than a few tens of minutes.

write_params() {
  id="$1"
  linear=false; longitudinal=false; survival=false
  covar_numeric="SEX age_at_baseline PC1 PC2 PC3"
  covar_categorical="site specimen"
  skip_pop_split=false
  mh_plot=false
  time_col="study_days"
  covarfile='${projectDir}/example/covariates.tsv'

  case "$id" in
    vcf_glm)
      input='${projectDir}/example/genotype/chr2[0-2].vcf'
      cache_id="smoke_vcf"; linear=true
      phenofile='${projectDir}/example/phenotype.cs2.tsv'; pheno_name="y,y2"
      interact=""; mh_plot=true ;;
    vcf_cph_interact)
      input='${projectDir}/example/genotype/chr2[0-2].vcf'
      cache_id="smoke_vcf"; survival=true
      phenofile='${projectDir}/example/phenotype.surv.tsv'; pheno_name="surv_y"
      interact="age_at_baseline"
      covarfile="$FIXTURES/covariates.nofid.tsv" ;;
    vcf_gallop)
      input='${projectDir}/example/genotype/chr2[0-2].vcf'
      cache_id="smoke_vcf"; longitudinal=true
      phenofile='${projectDir}/example/phenotype.lt.tsv'; pheno_name="y"
      # Set on purpose: GALLOP has no interaction support, so this asserts the
      # parameter is ignored rather than partially applied.
      interact="age_at_baseline" ;;
    plink_fid_glm_interact)
      input="$FIXTURES/plink_fid/chr2[0-2].pgen"
      cache_id="smoke_plink_fid"; linear=true
      phenofile='${projectDir}/example/phenotype.cs2.tsv'; pheno_name="y,y2"
      interact="age_at_baseline"
      covarfile="$FIXTURES/covariates.nofid.tsv" ;;
    plink_cph)
      input='${projectDir}/example/genotype_plink/chr2[0-2].pgen'
      cache_id="smoke_plink"; survival=true
      phenofile='${projectDir}/example/phenotype.surv.tsv'; pheno_name="surv_y"
      interact="" ;;
    plink_gallop)
      input='${projectDir}/example/genotype_plink/chr2[0-2].pgen'
      cache_id="smoke_plink"; longitudinal=true
      phenofile='${projectDir}/example/phenotype.lt.tsv'; pheno_name="y"
      interact=""
      covarfile="$FIXTURES/covariates.nofid.tsv" ;;
    bed_glm)
      input="$FIXTURES/plink_bed/chr2[0-2].bed"
      cache_id="smoke_bed"; linear=true
      phenofile='${projectDir}/example/phenotype.cs2.tsv'; pheno_name="y,y2"
      interact="" ;;
    skip_pop_split_glm)
      input='${projectDir}/example/genotype/chr2[0-2].vcf'
      cache_id="smoke_vcf_skip"; linear=true
      phenofile='${projectDir}/example/phenotype.cs2.tsv'; pheno_name="y,y2"
      interact=""; skip_pop_split=true ;;
    *) fail_hard "unknown scenario '$id'" ;;
  esac

  # covar_interact must also appear in covar_numeric (main.nf validates this).
  cat > "$PARAMSDIR/$id.yml" <<EOF
# Generated by tests/smoke_test.sh -- scenario: $id
# $(scenario_desc "$id")
input: "$input"
covarfile: "$covarfile"
phenofile: "$phenofile"

pheno_name: "$pheno_name"
covar_numeric: "$covar_numeric"
covar_categorical: "$covar_categorical"
covar_interact: "$interact"
covar_cat_min_count: 20
study_arm_col: "study_arm"
time_col: "$time_col"

linear_flag: $linear
longitudinal_flag: $longitudinal
survival_flag: $survival

skip_pop_split: $skip_pop_split
ancestry: "EUR"
assembly: "hg19"

r2thres: -9
geno: "0.05"
minor_allele_freq: "0.05"
minor_allele_ct: "20"
hwe: "1e-6"
mind: "0.05"
kinship: "0.177"

chunk_size: 30000
mh_plot: $mh_plot
publish_gwas_results: true

STORE_ROOT: "$STORE_ROOT"
PROJECT_NAME: "$PROJECT_NAME"
reference_dir: "$REFERENCE_DIR"
genetic_data_id: "$cache_id"
analysis_name: "$id"
EOF

  # Carried into the assertion step.
  SC_CACHE_ID="$cache_id"
  SC_INTERACT="$interact"
  SC_PHENOS="$(echo "$pheno_name" | tr ',' ' ')"
  if [ "$longitudinal" = true ]; then
    SC_MODEL="lmm_gallop"
  elif [ "$survival" = true ]; then
    SC_MODEL="cph"
  else
    SC_MODEL="glm"
  fi
}

# -------------------------------------------------------------- assertions
PASS=0
FAIL=0
REPORT_FILE=""

note_pass() { echo "PASS  $1  ($2)" >> "$REPORT_FILE"; PASS=$((PASS+1)); }
note_fail() {
  echo "FAIL  $1" >> "$REPORT_FILE"
  echo "$2" | sed '/^[[:space:]]*$/d; s/^/        - /' >> "$REPORT_FILE"
  FAIL=$((FAIL+1))
}

# 1-based index of a named column in a TSV header, or empty.
col_idx() { awk -F'\t' -v want="$2" 'NR==1{for(i=1;i<=NF;i++) if($i==want){print i; exit}}' "$1"; }

# Succeeds if the named column has at least one row that is not NA/empty.
has_value() {
  _idx="$(col_idx "$1" "$2")"
  [ -n "$_idx" ] || return 1
  awk -F'\t' -v i="$_idx" \
    'NR>1 && $i!="NA" && $i!="" && $i!="nan" {found=1; exit} END{exit !found}' "$1"
}

# Succeeds if any of the named columns is present.
any_col() {
  _f="$1"; shift
  for _c in "$@"; do
    [ -n "$(col_idx "$_f" "$_c")" ] && return 0
  done
  return 1
}

# Print any column name that appears more than once in the header.
dup_cols() {
  head -1 "$1" | tr '\t' '\n' | sort | uniq -d | tr '\n' ' '
}

NL='
'
check_scenario() {
  id="$1"
  errs=""
  add_err() { errs="$errs$1$NL"; }

  resdir="$PROJECT_DIR/analyses/$SC_CACHE_ID/$id/gwas_results/$SC_MODEL"
  if [ ! -d "$resdir" ]; then
    note_fail "$id" "no results directory at $resdir"
    return
  fi

  allres="$(find "$resdir" -maxdepth 1 -name '*_allresults.tsv' | sort)"
  if [ -z "$allres" ]; then
    note_fail "$id" "no *_allresults.tsv under $resdir"
    return
  fi

  # Every requested phenotype must have produced at least one result file.
  for p in $SC_PHENOS; do
    echo "$allres" | grep -q "_${p}_allresults\.tsv$" ||
      add_err "no results for phenotype '$p'"
  done

  for f in $allres; do
    base="$(basename "$f")"
    if [ "$(wc -l < "$f")" -lt 2 ]; then
      add_err "$base: header only, no result rows"
      continue
    fi

    # Harmonised schema, shared by GLM / GALLOP / CoxPH output.
    for c in '#CHROM' POS ID A1 BETA SE P P_INT P_2DF INTERACTION MODEL; do
      [ -n "$(col_idx "$f" "$c")" ] || add_err "$base: missing column '$c'"
    done

    has_value "$f" P || add_err "$base: every P value is NA"

    dups="$(dup_cols "$f")"
    [ -z "$dups" ] || add_err "$base: duplicated output columns: $dups"

    if [ "$SC_MODEL" = "lmm_gallop" ]; then
      # GALLOP reuses the *_INT columns for the SNP x TIME slope term -- that is
      # the point of the model, and it is unrelated to covar_interact. So they
      # must always be populated and INTERACTION must read "TIME", whether or
      # not covar_interact was set. Scenario vcf_gallop sets covar_interact on
      # purpose; this is where "GALLOP ignores it" gets asserted.
      has_value "$f" P_INT || add_err "$base: GALLOP produced no SNPxTIME term (P_INT all NA)"
      iidx="$(col_idx "$f" INTERACTION)"
      awk -F'\t' -v i="$iidx" 'NR>1 && $i!="TIME" {bad=1; exit} END{exit bad?1:0}' "$f" ||
        add_err "$base: GALLOP INTERACTION column is not 'TIME' (covar_interact leaked in?)"
    elif [ -n "$SC_INTERACT" ]; then
      has_value "$f" P_INT || add_err "$base: interaction requested but every P_INT is NA"
      iidx="$(col_idx "$f" INTERACTION)"
      awk -F'\t' -v i="$iidx" -v want="$SC_INTERACT" \
        'NR>1 && $i!=want {bad=1; exit} END{exit bad?1:0}' "$f" ||
        add_err "$base: INTERACTION column is not '$SC_INTERACT'"
    else
      # No interaction requested: the columns must exist but stay empty.
      if has_value "$f" P_INT; then
        add_err "$base: P_INT is populated but no interaction should have been applied"
      fi
    fi
  done

  # GLM only: one run must yield both a linear and a logistic fit, because y is
  # continuous and y2 is binary 0/1. plink2 names the statistic column
  # T_STAT/Z_STAT normally and T_OR_F_STAT/Z_OR_F_STAT under --glm interaction.
  if [ "$SC_MODEL" = "glm" ]; then
    case " $SC_PHENOS " in
      *" y2 "*)
        lin="$(find "$resdir/split" -name '*.y.results' 2>/dev/null | head -1)"
        log="$(find "$resdir/split" -name '*.y2.results' 2>/dev/null | head -1)"
        { [ -n "$lin" ] && any_col "$lin" T_STAT T_OR_F_STAT; } ||
          add_err "phenotype y did not produce a linear fit (no T_STAT/T_OR_F_STAT column)"
        { [ -n "$log" ] && any_col "$log" Z_STAT Z_OR_F_STAT; } ||
          add_err "phenotype y2 did not produce a logistic fit (no Z_STAT/Z_OR_F_STAT column)"
        ;;
    esac
  fi

  # The shared data-prep stage must have produced a harmonised covariate table.
  prep="$PROJECT_DIR/analyses/$SC_CACHE_ID/$id/prepared_data"
  if [ -z "$(find "$prep" -maxdepth 1 -name '*_filtered.pca.harmonized.tsv' 2>/dev/null)" ]; then
    add_err "no *_filtered.pca.harmonized.tsv in prepared_data"
  fi

  n_res="$(echo "$allres" | wc -l | tr -d ' ')"
  if [ -n "$errs" ]; then
    note_fail "$id" "$errs"
  else
    note_pass "$id" "$n_res result file(s), model=$SC_MODEL"
  fi
}

# ------------------------------------------------------------------- runner
run_scenario() {
  id="$1"
  write_params "$id"

  set -- nextflow -log "$LOGDIR/$id.nextflow.log" \
         run main.nf -profile "$PROFILE" \
         -params-file "$PARAMSDIR/$id.yml" \
         -work-dir "$PROJECT_DIR/work" \
         -ansi-log false
  [ -n "$CONTAINER" ] && set -- "$@" -process.container "$CONTAINER"

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "$*"
    return 0
  fi

  echo "==> [$id] $(scenario_desc "$id")"
  start=$SECONDS
  if "$@" > "$LOGDIR/$id.out" 2>&1; then
    echo "    nextflow ok in $((SECONDS-start))s"
    check_scenario "$id"
  else
    echo "    nextflow FAILED in $((SECONDS-start))s"
    tail -30 "$LOGDIR/$id.out" | sed 's/^/    | /'
    note_fail "$id" "nextflow exited non-zero (see $LOGDIR/$id.out)"
  fi
}

build_fixtures
REPORT_FILE="$LOGDIR/report.txt"
: > "$REPORT_FILE"

[ "$DRY_RUN" -eq 1 ] || echo "==> output directory: $OUTDIR"

SUITE_START=$SECONDS
for id in $SCENARIOS; do
  run_scenario "$id"
done

[ "$DRY_RUN" -eq 1 ] && exit 0

# -------------------------------------------------------------------- report
echo
echo "================= smoke test summary ================="
cat "$REPORT_FILE"
echo "------------------------------------------------------"
echo "passed: $PASS   failed: $FAIL   elapsed: $((SECONDS-SUITE_START))s"
echo "artifacts: $OUTDIR"
echo "======================================================"

if [ "$FAIL" -eq 0 ] && [ "$KEEP" -eq 0 ]; then
  echo "(removing $OUTDIR -- pass --keep to retain it)"
  rm -rf "$OUTDIR"
fi

[ "$FAIL" -eq 0 ] || exit 1
exit 0
