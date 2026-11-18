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
 *   nextflow run run_power_analysis.nf -params-file example/params_YML/power_before.yml
 *
 * Three modes:
 *   planning  BEFORE a run: --power_n_grid "500 1000 5000", no data files at all
 *   data      BEFORE a run: --phenofile + --covarfile, using the cohort you have
 *   after     AFTER a run: --power_analysis_set pointing at that run's
 *             prepared_data/*_filtered.pca.harmonized.tsv. This is the accurate one --
 *             N is the post-QC analysed sample count and the covariates are the real
 *             ones, PC1..PCn included, so the SD is a genuine residual SD.
 *
 * Required (data mode):
 *   --phenofile   Phenotype file (IID or #IID column)
 *   --covarfile   Covariate file, one row per subject (omit if --power_analysis_set is set)
 *   --pheno_name  Phenotype column(s); comma- or space-separated for several
 *   plus at least one analysis flag: --linear_flag (GLM), --longitudinal_flag (LMM)
 *   or --survival_flag (CoxPH).
 *
 * Optional (power-specific):
 *   --power_outdir       Output directory (default: "${launchDir}/power_analysis")
 *   --power_prefix       Output prefix (default: params.analysis_name, else "POWER")
 *   --power_alpha        Significance threshold (default: 5e-8)
 *   --power_target       Target power (default: 0.80)
 *   --power_maf_grid     MAF grid for the MDE columns (default: "0.01 0.05 0.10 0.20 0.50")
 *   --power_info         Imputation r2 applied to Var(G) (default: 1.0, i.e. assumed)
 *   --power_analysis_set Post-QC analysis set(s) from a completed run; a path, a glob or a
 *                        comma-separated list. Replaces --covarfile, so N, the covariates
 *                        and the PCs are the ones the GWAS actually used
 *   --power_n_grid       Planning mode: sample sizes to evaluate, e.g. "500 1000 5000"
 *   --power_icc          LMM: intraclass correlation; estimated from data if omitted,
 *                        required in planning mode
 *   --power_visits       LMM planning: mean measurements per subject (k)
 *   --power_event_rate   CoxPH planning: fraction with an event, so D = N * rate.
 *                        Without it, --power_n_grid is read as event counts
 *   --power_event_col    CoxPH: 0/1 event indicator column; defaults to pheno_name
 *   --power_slope_term   LMM: also report the SNPxtime slope bound (default: true)
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
params.power_analysis_set = params.containsKey('power_analysis_set') ? params.power_analysis_set : null
params.power_n_grid   = params.containsKey('power_n_grid')   ? params.power_n_grid   : null
params.power_icc        = params.containsKey('power_icc')        ? params.power_icc        : null
params.power_visits     = params.containsKey('power_visits')     ? params.power_visits     : null
params.power_event_rate = params.containsKey('power_event_rate') ? params.power_event_rate : null
params.power_event_col  = params.containsKey('power_event_col')  ? params.power_event_col  : null
params.power_slope_term = params.containsKey('power_slope_term') ? params.power_slope_term : true
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

// Cache per (path, delim) so multiple phenotype columns do not re-read the file.
// Value is Map<String, Integer> of { columnName -> nonBinaryCount }.
def __nonBinaryCountsCache = [:]

def countNonBinaryValues(path, colName, delim) {
    def cacheKey = "${path}||${delim}"
    def cached = __nonBinaryCountsCache[cacheKey]
    if (cached != null) {
        return (cached[colName] ?: 0) as int
    }

    def validValues = ['0', '1', '0.0', '1.0', '', 'NA', 'NaN'] as Set
    def countsByCol = [:]

    file(path).withReader { reader ->
        def headerLine = reader.readLine()
        if (!headerLine) {
            __nonBinaryCountsCache[cacheKey] = countsByCol
            return
        }

        def header = headerLine.split(delim).collect { it.trim().replaceFirst(/^#/, '') }
        // Initialize counts to 0 for each column
        header.each { h -> countsByCol[h] = 0 }

        String line
        while ((line = reader.readLine()) != null) {
            if (!line.trim()) continue
            def fields = line.split(delim, -1)

            // Count non-binary values for all columns in one pass
            for (int i = 0; i < header.size() && i < fields.size(); i++) {
                def v = fields[i].trim()
                if (!(v in validValues)) {
                    def key = header[i]
                    countsByCol[key] = (countsByCol[key] as int) + 1
                }
            }
        }
    }

    __nonBinaryCountsCache[cacheKey] = countsByCol
    return (countsByCol[colName] ?: 0) as int
}

/*
 * Validate up front so failures surface before the DAG starts.
 */
def nGrid = params.power_n_grid == null ? '' :
    ((params.power_n_grid instanceof List) ? params.power_n_grid.join(' ')
      : params.power_n_grid.toString().replaceAll(/[\[\],]/, ' ').trim().replaceAll(/\s+/, ' '))
def planningMode = nGrid as Boolean

// A glob is the natural way to name every study arm a run produced, so expand one
// here: file() returns a List for a pattern and a single Path otherwise.
// One YAML serves both modes, so refuse the ambiguous case rather than silently
// preferring one: the user must comment out the block they are not using.
if (planningMode && params.power_analysis_set) {
    error("--power_n_grid (planning, BEFORE a run) and --power_analysis_set (AFTER a run) " +
          "are mutually exclusive.\nComment out whichever block you are not using.")
}

def analysisSetFiles = []
if (params.power_analysis_set && !planningMode) {
    params.power_analysis_set.toString().split(',').collect { it.trim() }.findAll { it }.each { pat ->
        def matched = file(pat)
        if (matched instanceof List) analysisSetFiles.addAll(matched)
        else if (matched.exists())   analysisSetFiles << matched
        else error("--power_analysis_set matched no file: '${pat}'")
    }
    if (!analysisSetFiles) {
        error("--power_analysis_set '${params.power_analysis_set}' matched no files.\n" +
              "Expected a completed run's prepared_data/*_filtered.pca.harmonized.tsv.")
    }
    analysisSetFiles = analysisSetFiles.sort { it.name }
}
def useAnalysisSet = analysisSetFiles as Boolean

def missingParams = []
if (!planningMode) {
    if (!params.phenofile) missingParams << '--phenofile (path to phenotype file)'
    if (!params.covarfile && !useAnalysisSet) {
        missingParams << '--covarfile (path to covariates file), or --power_analysis_set for a completed run'
    }
}
if (missingParams) {
    error("Missing required parameter(s):\n" +
          missingParams.collect { "  - ${it}" }.join('\n') +
          "\nGenotype files are NOT required for power analysis." +
          "\nFor study planning without any data, set --power_n_grid \"500 1000 5000\" instead." +
          "\nSee example/params_YML/power_before.yml for a template.")
}

def designs = []
if (params.linear_flag)       designs << 'glm'
if (params.longitudinal_flag) designs << 'lmm'
if (params.survival_flag)     designs << 'cph'
if (!designs) {
    error("Set at least one analysis flag: --linear_flag (GLM), --longitudinal_flag (LMM) " +
          "or --survival_flag (CoxPH).")
}
def phenoNames = planningMode ? ['planned'] : parsePhenoNames(params.pheno_name)
def phenoDelim  = planningMode ? '\t' : detectDelimiter(params.phenofile)
def phenoHeader = planningMode ? ([] as Set) : readHeaderColumns(params.phenofile, '--phenofile')
def covarLabel  = useAnalysisSet ? '--power_analysis_set' : '--covarfile'
// Keep the Path object: toString() drops the gs:// scheme on cloud paths.
def covarSource = useAnalysisSet ? analysisSetFiles[0] : params.covarfile
def covarHeader = planningMode ? ([] as Set) : readHeaderColumns(covarSource, covarLabel)

def missingColumns = []
if (!planningMode) {
if (!phenoNames) missingColumns << "--pheno_name is empty"
phenoNames.each { col ->
    if (!(col in phenoHeader)) missingColumns << "--pheno_name '${col}' not found in --phenofile columns"
}
// PCs are computed later by COMPUTE_PCA, so a raw covariate file legitimately lacks
// them -- but a post-QC analysis set must have them, so check for them there.
def numericList = (params.covar_numeric ?: '').split(/[\s,]+/).findAll { it }
numericList.each { col ->
    def pcExempt = !useAnalysisSet && (col ==~ /(?i)^PC\d+$/)
    if (!pcExempt && !(col in covarHeader)) {
        missingColumns << "--covar_numeric '${col}' not found in ${covarLabel} columns"
    }
}
(params.covar_categorical ?: '').split(/[\s,]+/).findAll { it }.each { col ->
    if (!(col in covarHeader)) missingColumns << "--covar_categorical '${col}' not found in ${covarLabel} columns"
}
if (params.study_arm_col && !(params.study_arm_col in covarHeader)) {
    missingColumns << "--study_arm_col '${params.study_arm_col}' not found in ${covarLabel} columns"
}
if ('lmm' in designs && !(params.time_col in phenoHeader)) {
    missingColumns << "--time_col '${params.time_col}' not found in --phenofile columns (required for LMM)"
}
// GLM/LMM need a quantitative phenotype; CoxPH needs a 0/1 event indicator.
phenoNames.each { col ->
    if (col in phenoHeader) {
        def isBinary = countNonBinaryValues(params.phenofile, col, phenoDelim) == 0
        if (isBinary && (('glm' in designs) || ('lmm' in designs))) {
            missingColumns << ("--pheno_name '${col}' looks binary (only 0/1 values); binary/logistic " +
                               "power is not supported, since N_info for a case/control contrast is " +
                               "n1*n2/(n1+n2), not N")
        }
        if (!isBinary && 'cph' in designs && !params.power_event_col) {
            missingColumns << ("--pheno_name '${col}' is not 0/1, but CoxPH needs an event indicator; " +
                               "set --power_event_col to the 0/1 column")
        }
    }
}
}
if (missingColumns) {
    error("Phenotype/covariate validation failed:\n" +
          missingColumns.collect { "  - ${it}" }.join('\n') +
          "\n\nphenofile ('${params.phenofile}') columns found: ${phenoHeader.sort()}" +
          "\n${covarLabel} ('${covarSource}') columns found: ${covarHeader.sort()}")
}

// A YAML list arrives as a Groovy List, which would interpolate as "[0.01, 0.05]".
def mafGrid = (params.power_maf_grid instanceof List)
    ? params.power_maf_grid.join(' ')
    : params.power_maf_grid.toString().replaceAll(/[\[\],]/, ' ').trim().replaceAll(/\s+/, ' ')

log.info """\
 PROGWAS - POWER ANALYSIS
 ======================================
 mode                 : ${planningMode ? 'planning (--power_n_grid)' : (useAnalysisSet ? 'after a run (post-QC analysis set)' : 'data (phenotype + covariate files)')}
 phenotype file       : ${planningMode ? '(none)' : params.phenofile}
 covariate file       : ${planningMode ? '(none)' : (useAnalysisSet ? '(superseded by the analysis set)' : params.covarfile)}
 analysis set         : ${planningMode ? '(not applicable in planning mode)' : (useAnalysisSet ? analysisSetFiles.collect { it.name }.join(', ') : '(none - N is PRE-QC and therefore optimistic)')}
 phenotype name(s)    : ${planningMode ? "N grid: ${nGrid}" : phenoNames.join(', ')}
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
 * One process per design. They differ only in the script called and the
 * design-specific flags; N_info is what actually separates them.
 */
// Shell-quote a value for the task command, escaping embedded single quotes.
// Staged path inputs (phenofile, covarfile) are already escaped by Nextflow, so they are not wrapped again.
def q = { v -> "'" + v.toString().replace("'", "'\\''") + "'" }

def commonOpts = { pheno ->
    """--pheno-name ${q(pheno)} \\
        --study-arm-col ${q(params.study_arm_col ?: '')} \\
        --covar-numeric ${q(params.covar_numeric ?: '')} \\
        --covar-categorical ${q(params.covar_categorical ?: '')}"""
}
// The analysis set supersedes the covariate file when present. Staged names are
// used, not the original paths, since the task runs in its own work directory.
def covarOpts = { covarfile, asetfiles ->
    def names = [asetfiles].flatten().findAll { it }.collect { it.name }
    names ? "--analysis-set ${q(names.join(','))}" : "--covar-file ${covarfile}"
}

def sharedOpts = {
    def sd = params.power_sd_pheno ? "--sd-pheno ${params.power_sd_pheno}" : ''
    def rank = params.power_rank_int ? '--rank-int' : ''
    """--maf-grid ${q(mafGrid)} \\
        --minor-allele-freq ${params.minor_allele_freq} \\
        --minor-allele-ct ${params.minor_allele_ct} \\
        --alpha ${params.power_alpha} \\
        --target-power ${params.power_target} \\
        --info ${params.power_info} \\
        --residualize ${params.power_residualize} \\
        --per-arm ${params.power_per_arm} \\
        ${sd} ${rank}"""
}

process POWER_GLM {
    tag "${pheno}"
    label 'small'
    errorStrategy 'terminate'

    publishDir "${params.power_outdir}/glm/${params.power_prefix}", mode: 'copy',
               overwrite: true, pattern: "*_power_{results,extended_info}.tsv"
    publishDir "${params.power_outdir}/glm/${params.power_prefix}", mode: 'copy',
               overwrite: true, pattern: "{power_results,extended_info}_key.tsv"
    publishDir "${params.power_outdir}/glm/${params.power_prefix}", mode: 'copy',
               overwrite: true, pattern: "*_power_curve.png"
    publishDir "${params.power_outdir}/glm/${params.power_prefix}/provenance", mode: 'copy',
               overwrite: true, pattern: "*_provenance.txt"

    input:
      tuple val(pheno), path(phenofile), path(covarfile), path(analysisset)

    output:
      path "${params.power_prefix}_${pheno}_power_results.tsv",       emit: results
      path "${params.power_prefix}_${pheno}_power_extended_info.tsv", emit: extended
      path "power_results_key.tsv",                                   emit: results_key
      path "extended_info_key.tsv",                                   emit: extended_key
      path "${params.power_prefix}_${pheno}_power_curve.png",         emit: plot, optional: true
      path "${params.power_prefix}_${pheno}_provenance.txt",          emit: provenance

    script:
      def dataOpts = planningMode ? "--n-grid ${q(nGrid)}"
                                  : "--pheno-file ${phenofile} \\\n        ${covarOpts(covarfile, analysisset)} \\\n        ${commonOpts(pheno)}"
      """
      power_analysis_glm.R \\
        ${dataOpts} \\
        ${sharedOpts()} \\
        --out-prefix ${q("${params.power_prefix}_${pheno}")}
      """
}

process POWER_LMM {
    tag "${pheno}"
    label 'small'
    errorStrategy 'terminate'

    publishDir "${params.power_outdir}/lmm_gallop/${params.power_prefix}", mode: 'copy',
               overwrite: true, pattern: "*_power_{results,extended_info}.tsv"
    publishDir "${params.power_outdir}/lmm_gallop/${params.power_prefix}", mode: 'copy',
               overwrite: true, pattern: "{power_results,extended_info}_key.tsv"
    publishDir "${params.power_outdir}/lmm_gallop/${params.power_prefix}", mode: 'copy',
               overwrite: true, pattern: "*_power_curve.png"
    publishDir "${params.power_outdir}/lmm_gallop/${params.power_prefix}/provenance", mode: 'copy',
               overwrite: true, pattern: "*_provenance.txt"

    input:
      tuple val(pheno), path(phenofile), path(covarfile), path(analysisset)

    output:
      path "${params.power_prefix}_${pheno}_power_results.tsv",       emit: results
      path "${params.power_prefix}_${pheno}_power_extended_info.tsv", emit: extended
      path "power_results_key.tsv",                                   emit: results_key
      path "extended_info_key.tsv",                                   emit: extended_key
      path "${params.power_prefix}_${pheno}_power_curve.png",         emit: plot, optional: true
      path "${params.power_prefix}_${pheno}_provenance.txt",          emit: provenance

    script:
      def iccOpt = params.power_icc != null ? "--icc ${params.power_icc}" : ''
      def visitsOpt = params.power_visits != null ? "--visits ${params.power_visits}" : ''
      def dataOpts = planningMode
          ? "--n-grid ${q(nGrid)} ${visitsOpt}"
          : "--pheno-file ${phenofile} \\\n        ${covarOpts(covarfile, analysisset)} \\\n        --time-col ${q(params.time_col)} \\\n        ${commonOpts(pheno)}"
      """
      power_analysis_lmm.R \\
        ${dataOpts} \\
        ${iccOpt} \\
        --slope-term ${params.power_slope_term} \\
        ${sharedOpts()} \\
        --out-prefix ${q("${params.power_prefix}_${pheno}")}
      """
}

process POWER_CPH {
    tag "${pheno}"
    label 'small'
    errorStrategy 'terminate'

    publishDir "${params.power_outdir}/cph/${params.power_prefix}", mode: 'copy',
               overwrite: true, pattern: "*_power_{results,extended_info}.tsv"
    publishDir "${params.power_outdir}/cph/${params.power_prefix}", mode: 'copy',
               overwrite: true, pattern: "{power_results,extended_info}_key.tsv"
    publishDir "${params.power_outdir}/cph/${params.power_prefix}", mode: 'copy',
               overwrite: true, pattern: "*_power_curve.png"
    publishDir "${params.power_outdir}/cph/${params.power_prefix}/provenance", mode: 'copy',
               overwrite: true, pattern: "*_provenance.txt"

    input:
      tuple val(pheno), path(phenofile), path(covarfile), path(analysisset)

    output:
      path "${params.power_prefix}_${pheno}_power_results.tsv",       emit: results
      path "${params.power_prefix}_${pheno}_power_extended_info.tsv", emit: extended
      path "power_results_key.tsv",                                   emit: results_key
      path "extended_info_key.tsv",                                   emit: extended_key
      path "${params.power_prefix}_${pheno}_power_curve.png",         emit: plot, optional: true
      path "${params.power_prefix}_${pheno}_provenance.txt",          emit: provenance

    script:
      def rateOpt = params.power_event_rate != null ? "--event-rate ${params.power_event_rate}" : ''
      def evOpt = params.power_event_col ? "--event-col ${q(params.power_event_col)}" : ''
      def dataOpts = planningMode
          ? "--n-grid ${q(nGrid)} ${rateOpt}"
          : "--pheno-file ${phenofile} \\\n        ${covarOpts(covarfile, analysisset)} \\\n        ${evOpt} \\\n        ${commonOpts(pheno)}"
      """
      power_analysis_cph.R \\
        ${dataOpts} \\
        ${sharedOpts()} \\
        --out-prefix ${q("${params.power_prefix}_${pheno}")}
      """
}

workflow {
    // Planning mode stages no files; an empty list keeps the input cardinality the same.
    // The covariate file is likewise left unstaged when an analysis set supersedes it.
    pheno_f = planningMode ? [] : file(params.phenofile, checkIfExists: true)
    covar_f = (planningMode || useAnalysisSet) ? [] : file(params.covarfile, checkIfExists: true)
    aset_f  = analysisSetFiles

    jobs = Channel.fromList(phenoNames).map { p -> tuple(p, pheno_f, covar_f, aset_f) }

    // All requested designs run in one launch, unlike main.nf which picks exactly one:
    // power is cheap, and the designs are directly comparable once N_info is right.
    if ('glm' in designs) { POWER_GLM(jobs); POWER_GLM.out.results.subscribe { f -> log.info " glm:  ${f.name}" } }
    if ('lmm' in designs) { POWER_LMM(jobs); POWER_LMM.out.results.subscribe { f -> log.info " lmm:  ${f.name}" } }
    if ('cph' in designs) { POWER_CPH(jobs); POWER_CPH.out.results.subscribe { f -> log.info " cph:  ${f.name}" } }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
