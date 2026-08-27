#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Minimal METAL workflow for GWAS summary statistics.
 *
 * Required:
 *   --metal_input    Glob pattern or comma-separated list of input files
 *   long gwas output with the folloiwng columns (header names are case-sensitive):
 *  #CHROM	POS	ID	REF	ALT	A1	A1_FREQ	MISS_FREQ	OBS_CT	TEST	BETA	exp(BETA)	SE	P
 *  This pipeline adds A2 column based on A1 and REF/ALT
 *
 * Optional:
 *   --metal_outdir      Output directory (default: "${launchDir}/metal_results")
 *   --metal_prefix      Output prefix (default: "SURV_META")
 *   --metal_pheno_name  Comma-separated phenotype names. If set, input files are
 *                       grouped by phenotype suffix (_<phenotype>_allresults.tsv)
 *                       and one METAL run is launched per phenotype.
 *   --metal_effect      Which effect columns to meta-analyse (default: unset -> GLM/CPH):
 *                         unset       BETA / SE / P            (GWASGLM, GWASCPH)
 *                         slope       BETAs / SEs / Ps          (GWASGALLOP slope effect)
 *                         intercept   BETAi / SEi / Pi          (GWASGALLOP intercept effect)
 */

params.metal_input      = params.metal_input ?: null
params.metal_outdir     = params.metal_outdir ?: "${launchDir}/metal_results"
params.metal_prefix     = params.metal_prefix ?: "SURV_META"
params.metal_pheno_name = params.metal_pheno_name ?: null
params.metal_effect     = params.metal_effect ?: null

if (!params.metal_input) {
    error "Missing required parameter: --metal_input"
}

def metalEffectColumns = [
    'slope'    : ['BETAs', 'SEs', 'Ps'],
    'intercept': ['BETAi', 'SEi', 'Pi'],
].withDefault { null }

if (params.metal_effect && !metalEffectColumns.containsKey(params.metal_effect)) {
    error "Invalid --metal_effect '${params.metal_effect}': expected 'slope' or 'intercept' (omit for GLM/CPH default, which uses BETA/SE/P)"
}

def (metalEffectCol, metalStderrCol, metalPvalCol) = metalEffectColumns[params.metal_effect] ?: ['BETA', 'SE', 'P']

def parsePhenoNames(value) {
  (value ?: '').split(',').collect { it.trim() }.findAll { it }
}

def matchPhenoName(filename, names) {
  def base = filename.replaceAll(/\.gz$/, '')
  def matches = names.findAll { base.endsWith("_${it}_allresults.tsv") }
  matches ? matches.max { it.length() } : null
}

def phenoNames = parsePhenoNames(params.metal_pheno_name)

def inputPatterns = params.metal_input
    .toString()
    .split(',')
    .collect { it.trim() }
    .findAll { it }

Channel
    .fromPath(inputPatterns, checkIfExists: true)
    .ifEmpty { error "No files matched --metal_input: ${params.metal_input}" }
  .map { f ->
    def pheno = phenoNames ? matchPhenoName(f.getName(), phenoNames) : ''
    if (phenoNames && pheno == null) {
      error "File does not match any --metal_pheno_name value (${phenoNames.join(', ')}): ${f}"
    }
    tuple(pheno, f)
  }
  .groupTuple()
  .set { metal_input_groups_ch }

process RUN_METAL {
    label 'medium'
    publishDir "${params.metal_outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(pheno), path(sumstats_files, stageAs: 'input??/*')

    output:
    path "${prefix}_metal_script.txt", emit: script
    path "prepared/*.metal.tsv.gz", emit: prepared
    path "${prefix}_command.log", emit: log, optional: true
    path "${prefix}*", emit: results

    script:
    prefix = pheno ? "${params.metal_prefix}_${pheno}" : params.metal_prefix
    """
    set -euo pipefail

    mkdir -p prepared

    cat > ${prefix}_metal_script.txt << 'METAL_EOF'
SCHEME STDERR
GENOMICCONTROL OFF
AVERAGEFREQ ON
MINMAXFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize AS OBS_CT
MARKER ID
FREQ A1_FREQ
ALLELE A1 A2
EFFECT ${metalEffectCol}
STDERR ${metalStderrCol}
PVAL ${metalPvalCol}
WEIGHT OBS_CT
METAL_EOF

    total_ref=0
    total_alt=0

    for f in ${sumstats_files}; do
      base=\$(basename "\$f")
      tag=\$(basename "\$(dirname "\$f")")
      out="prepared/\${tag}_\${base%.gz}.metal.tsv"

      ( [[ "\$f" == *.gz ]] && zcat "\$f" || cat "\$f" ) | awk -v count_file="\${out}.counts" 'BEGIN{FS=OFS="\\t"; c_ref=0; c_alt=0}
        NR==1 {
          for (i=1; i<=NF; i++) {
            if (\$i == "REF") ref_col=i
            else if (\$i == "ALT") alt_col=i
            else if (\$i == "A1") a1_col=i
          }
          print \$0, "A2";
          next
        }
        {
          a2="."
          if (\$a1_col == \$ref_col) { a2=\$alt_col; c_ref++ }
          else if (\$a1_col == \$alt_col) { a2=\$ref_col; c_alt++ }
          print \$0, a2
        }
        END { printf("%d\\t%d\\n", c_ref, c_alt) > count_file }' > "\$out"

      a1_eq_ref=\$(cut -f1 "\${out}.counts")
      a1_eq_alt=\$(cut -f2 "\${out}.counts")
      rm -f "\${out}.counts"

      total_ref=\$((total_ref + a1_eq_ref))
      total_alt=\$((total_alt + a1_eq_alt))

      echo "[A2_BUILD] file=\$base A1_EQ_REF=\$a1_eq_ref A1_EQ_ALT=\$a1_eq_alt"
      gzip -f "\$out"
      echo "PROCESS \${out}.gz" >> ${prefix}_metal_script.txt
    done
    
    echo "[A2_BUILD] total A1_EQ_REF=\$total_ref A1_EQ_ALT=\$total_alt"

    cat >> ${prefix}_metal_script.txt << METAL_EOF
  OUTFILE ${prefix} .TBL
ANALYZE HETEROGENEITY
QUIT
METAL_EOF

    metal ${prefix}_metal_script.txt
    
    cp .command.log ${prefix}_command.log || true
    """
}

workflow {
    RUN_METAL(metal_input_groups_ch)
}
