#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Study power analysis for proGWAS.
 *
 * Computes minimum detectable effect sizes for the analysis you are about to run
 * (or already ran), from the phenotype and covariate files. No genotypes, no
 * References and no plink2 are required, so this is cheap to run at planning time.
 *
 * Point it at an existing params file -- the data parameters are the same ones
 * main.nf uses:
 *   nextflow run run_power_analysis.nf -params-file example/params_YML/power_cs.yml
 *
 * Required:
 *   --phenofile   Phenotype file (IID or #IID column)
 *   --covarfile   Covariate file, one row per subject
 *   --pheno_name  Phenotype column(s); comma- or space-separated for several
 *   plus one analysis flag: --linear_flag (GLM). Longitudinal (LMM) and survival
 *   (CoxPH) power are not implemented yet.
 *
 * Optional (power-specific):
 *   --power_outdir       Output directory (default: "${launchDir}/power_analysis")
 *   --power_prefix       Output prefix (default: params.analysis_name, else "POWER")
 *   --power_alpha        Significance threshold (default: 5e-8)
 *   --power_target       Target power (default: 0.80)
 *   --power_maf_grid     MAF grid for the MDE columns (default: "0.01 0.05 0.10 0.20 0.50")
 *   --power_info         Imputation r2 applied to Var(G) (default: 1.0)
 *   --power_sd_pheno     Measured phenotype SD; overrides the computed value
 *   --power_rank_int     Phenotype was rank-INT transformed (default: false)
 *   --power_residualize  Take the SD after regressing out covariates (default: true)
 *   --power_per_arm      Report each study arm separately (default: true)
 *
 * Reused from the main pipeline: --covar_numeric, --covar_categorical,
 * --study_arm_col, --minor_allele_freq, --minor_allele_ct, --analysis_name.
 */

// containsKey rather than '?:' throughout: '?:' turns an explicit false or 0 back
// into the default, and reading an unset param emits an "undefined parameter" warning.
params.power_outdir   = params.containsKey('power_outdir')   ? params.power_outdir   : "${launchDir}/power_analysis"
params.power_prefix   = params.containsKey('power_prefix')   ? params.power_prefix   : (params.analysis_name ?: 'POWER')
params.power_alpha    = params.containsKey('power_alpha')    ? params.power_alpha    : 5e-8
params.power_target   = params.containsKey('power_target')   ? params.power_target   : 0.80
params.power_maf_grid = params.containsKey('power_maf_grid') ? params.power_maf_grid : '0.01 0.05 0.10 0.20 0.50'
params.power_info     = params.containsKey('power_info')     ? params.power_info     : 1.0
params.power_sd_pheno = params.containsKey('power_sd_pheno') ? params.power_sd_pheno : null
params.power_rank_int    = params.containsKey('power_rank_int')    ? params.power_rank_int    : false
params.power_residualize = params.containsKey('power_residualize') ? params.power_residualize : true
params.power_per_arm     = params.containsKey('power_per_arm')     ? params.power_per_arm     : true

/*
 * Validation helpers, mirroring main.nf so the same params file behaves the same way.
 */
def readHeaderColumns(path, label) {
    def line
    try {
        line = file(path).withReader { it.readLine() }
    } catch (Exception e) {
        error("Could not read ${label} '${path}': ${e.message}")
    }
    if (!line) {
        error("${label} '${path}' appears to be empty (no header line found).")
    }
    def delim = line.contains('\t') ? '\t' : ','
    return line.split(delim).collect { it.trim().replaceFirst(/^#/, '') } as Set
}

def detectDelimiter(path) {
    def line = file(path).withReader { it.readLine() }
    return line?.contains('\t') ? '\t' : ','
}

def parsePhenoNames(value) {
    (value ?: '').split(/[\s,]+/).collect { it.trim() }.findAll { it }
}

def countNonBinaryValues(path, colName, delim) {
    def lines = file(path).readLines()
    if (lines.size() < 2) return 0
    def header = lines[0].split(delim).collect { it.trim().replaceFirst(/^#/, '') }
    def idx = header.findIndexOf { it == colName }
    if (idx < 0) return 0
    def validValues = ['0', '1', '0.0', '1.0', '', 'NA', 'NaN'] as Set
    return lines.drop(1).count { line ->
        if (!line.trim()) return false
        def fields = line.split(delim, -1)
        idx < fields.size() && !(fields[idx].trim() in validValues)
    }
}

/*
 * Validate up front so failures surface before the DAG starts.
 */
def missingParams = []
if (!params.phenofile) missingParams << '--phenofile (path to phenotype file)'
if (!params.covarfile) missingParams << '--covarfile (path to covariates file)'
if (missingParams) {
    error("Missing required parameter(s):\n" +
          missingParams.collect { "  - ${it}" }.join('\n') +
          "\nGenotype files are NOT required for power analysis." +
          "\nSee example/params_YML/power_cs.yml for a template.")
}

def designs = []
if (params.linear_flag)       designs << 'glm'
if (params.longitudinal_flag) designs << 'lmm'
if (params.survival_flag)     designs << 'cph'
if (!designs) {
    error("Set at least one analysis flag: --linear_flag (GLM). " +
          "Longitudinal (LMM) and survival (CoxPH) power are not implemented yet.")
}
def unsupported = designs - ['glm']
if (unsupported) {
    error("Power analysis is only implemented for GLM (--linear_flag) so far; requested " +
          "${unsupported.join(', ')}.\nN_info differs by design -- N for GLM, N*k/(1+(k-1)*rho) " +
          "for LMM, and events D for CoxPH -- so the GLM numbers must not be reused for them.")
}

def phenoDelim = detectDelimiter(params.phenofile)
def phenoHeader = readHeaderColumns(params.phenofile, '--phenofile')
def covarHeader = readHeaderColumns(params.covarfile, '--covarfile')
def phenoNames = parsePhenoNames(params.pheno_name)

def missingColumns = []
if (!phenoNames) missingColumns << "--pheno_name is empty"
phenoNames.each { col ->
    if (!(col in phenoHeader)) missingColumns << "--pheno_name '${col}' not found in --phenofile columns"
}
// PCs are computed later by COMPUTE_PCA, so they are not expected in a raw covariate file.
def numericList = (params.covar_numeric ?: '').split(/[\s,]+/).findAll { it }
numericList.each { col ->
    if (!(col ==~ /(?i)^PC\d+$/) && !(col in covarHeader)) {
        missingColumns << "--covar_numeric '${col}' not found in --covarfile columns"
    }
}
(params.covar_categorical ?: '').split(/[\s,]+/).findAll { it }.each { col ->
    if (!(col in covarHeader)) missingColumns << "--covar_categorical '${col}' not found in --covarfile columns"
}
if (params.study_arm_col && !(params.study_arm_col in covarHeader)) {
    missingColumns << "--study_arm_col '${params.study_arm_col}' not found in --covarfile columns"
}
// A binary phenotype needs N_info = n1*n2/(n1+n2), not N, so reject it here.
phenoNames.each { col ->
    if (col in phenoHeader && countNonBinaryValues(params.phenofile, col, phenoDelim) == 0) {
        missingColumns << ("--pheno_name '${col}' looks binary (only 0/1 values); binary/logistic " +
                           "power is not supported, since N_info for a case/control contrast is " +
                           "n1*n2/(n1+n2), not N")
    }
}
if (missingColumns) {
    error("Phenotype/covariate validation failed:\n" +
          missingColumns.collect { "  - ${it}" }.join('\n') +
          "\n\nphenofile ('${params.phenofile}') columns found: ${phenoHeader.sort()}" +
          "\ncovarfile ('${params.covarfile}') columns found: ${covarHeader.sort()}")
}

// A YAML list arrives as a Groovy List, which would interpolate as "[0.01, 0.05]".
def mafGrid = (params.power_maf_grid instanceof List)
    ? params.power_maf_grid.join(' ')
    : params.power_maf_grid.toString().replaceAll(/[\[\],]/, ' ').trim().replaceAll(/\s+/, ' ')

log.info """\
 PROGWAS - POWER ANALYSIS
 ======================================
 phenotype file       : ${params.phenofile}
 covariate file       : ${params.covarfile}
 phenotype name(s)    : ${phenoNames.join(', ')}
 numeric covariates   : ${params.covar_numeric ?: '(none)'}
 categorical covars   : ${params.covar_categorical ?: '(none)'}
 study arm column     : ${params.study_arm_col ?: '(none)'}
 designs              : ${designs.join(', ')}
 alpha / target power : ${params.power_alpha} / ${params.power_target}
 MAF grid             : ${mafGrid}
 MAF / MAC thresholds : ${params.minor_allele_freq} / ${params.minor_allele_ct}
 imputation r2        : ${params.power_info}${params.power_info == 1.0 ? ' (assumed; no genotypes read)' : ''}
 output directory     : ${params.power_outdir}
 """

/*
 * Process: GLM power analysis, one task per phenotype.
 */
process POWER_GLM {
    tag "${pheno}"
    label 'small'
    // conf/base.config retries twice by default, which would bury a deterministic argument error.
    errorStrategy 'terminate'

    publishDir "${params.power_outdir}", mode: 'copy', overwrite: true,
               pattern: "*_power_{table,long}.tsv"
    publishDir "${params.power_outdir}/provenance", mode: 'copy', overwrite: true,
               pattern: "*_provenance.txt"

    input:
      tuple val(pheno), path(phenofile), path(covarfile)

    output:
      path "${params.power_prefix}_${pheno}_power_table.tsv", emit: table
      path "${params.power_prefix}_${pheno}_power_long.tsv",  emit: tidy
      path "${params.power_prefix}_${pheno}_provenance.txt",  emit: provenance

    script:
      def sdOpt = params.power_sd_pheno ? "--sd-pheno ${params.power_sd_pheno}" : ''
      def rankOpt = params.power_rank_int ? '--rank-int' : ''
      """
      power_analysis_glm.R \\
        --pheno-file ${phenofile} \\
        --covar-file ${covarfile} \\
        --pheno-name '${pheno}' \\
        --study-arm-col '${params.study_arm_col ?: ''}' \\
        --covar-numeric '${params.covar_numeric ?: ''}' \\
        --covar-categorical '${params.covar_categorical ?: ''}' \\
        --maf-grid '${mafGrid}' \\
        --minor-allele-freq ${params.minor_allele_freq} \\
        --minor-allele-ct ${params.minor_allele_ct} \\
        --alpha ${params.power_alpha} \\
        --target-power ${params.power_target} \\
        --info ${params.power_info} \\
        --residualize ${params.power_residualize} \\
        --per-arm ${params.power_per_arm} \\
        ${sdOpt} ${rankOpt} \\
        --out-prefix ${params.power_prefix}_${pheno}
      """
}

workflow {
    pheno_f = file(params.phenofile, checkIfExists: true)
    covar_f = file(params.covarfile, checkIfExists: true)

    jobs = Channel.fromList(phenoNames).map { p -> tuple(p, pheno_f, covar_f) }

    POWER_GLM(jobs)

    POWER_GLM.out.table.subscribe { f -> log.info " power table: ${f.name}" }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
