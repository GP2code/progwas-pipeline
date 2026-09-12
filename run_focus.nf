#!/usr/bin/env nextflow

/*
 * SNP-Focused Analysis Workflow
 * Performs targeted GWAS on a subset of SNPs with optional stratification
 * 
 * Usage:
 *   nextflow run run_focus.nf -params-file example/params_YML/focus_cs.yml -profile standard
 */

nextflow.enable.dsl = 2

// Publish GWAS results when running this focused analysis workflow
params.publish_gwas_results = true

// ==================================================================================
// Log Configuration
// ==================================================================================
log.info """\
 SNP-FOCUSED GWAS ANALYSIS
 ==========================================
 Focus plink input                : ${params.focus_plink_input}
 Focus covariate file             : ${params.focus_covar_file}
 Focus phenotype name             : ${params.focus_pheno_name}
 Focus strata file                : ${params.focus_strata_file ?: 'None (single analysis)'}
 Phenotype file                   : ${params.phenofile}
 Longitudinal analysis            : ${params.longitudinal_flag}
 Survival analysis                : ${params.survival_flag}
 Numeric covariates               : ${params.covar_numeric}
 Categorical covariates           : ${params.covar_categorical}
 Interaction covariate            : ${params.covar_interact ?: 'None'}
 Analysis name                    : ${params.analysis_name}
 """

// ==================================================================================
// Import Processes from Main Modules
// ==================================================================================
include { HARMONIZE_CATEGORICAL_COVARS } from './modules/dataprep.nf'
include { EXPORT_PLINK; RAWFILE_EXPORT } from './modules/dataprep.nf'
include { GWASGLM; GWASGALLOP; GWASCPH } from './modules/gwas.nf'

// ==================================================================================
// Process: Rename covariate file to expected pattern
// ==================================================================================
process RENAME_COVAR_FILE {
  scratch true
  label 'small'
  
  input:
    path covar_file
  
  output:
    path "${params.ancestry}_focus_filtered.pca.tsv"
  
  script:
    """
    cp ${covar_file} ${params.ancestry}_focus_filtered.pca.tsv
    """
}

// ==================================================================================
// Process: Normalize variant IDs to chr:pos:ref:alt format
// Required for survival.R and gallop.py which parse IDs by splitting on ':'
// ==================================================================================
process NORMALIZE_VARIANT_IDS {
  scratch true
  label 'small'

  input:
    tuple val(fileTag), path(pgen), path(pvar), path(psam)

  output:
    tuple val("${fileTag}_normid"), path("${fileTag}_normid.pgen"), path("${fileTag}_normid.pvar"), path("${fileTag}_normid.psam")

  script:
    """
    plink2 \\
        --pfile ${fileTag} \\
        --set-all-var-ids @:#:\\\$r:\\\$a \\
        --new-id-max-allele-len 1000 \\
        --output-chr chrM \\
        --make-pgen \\
        --out ${fileTag}_normid \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toMega()}
    """
}

// ==================================================================================
// New Process: Filter Samples by STRATA
// ==================================================================================
process FILTER_BY_STRATA {
  scratch true
  label 'small'
  
  input:
    path samplelist
    path strata_file
    each strata
  
  output:
    tuple val(strata), path("*_filtered.pca.harmonized.tsv")
  
  script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    
    # Load sample list and strata file
    samples_df = pd.read_csv("${samplelist}", sep="\\t")
    strata_df = pd.read_csv("${strata_file}", sep="\\t")
    
    # Normalize FID column name if needed
    if '#FID' in strata_df.columns and 'FID' not in strata_df.columns:
        strata_df.rename(columns={'#FID': 'FID'}, inplace=True)
    
    # Filter strata file for this specific strata value
    strata_samples = strata_df[strata_df['STRATA'] == '${strata}']
    
    # Get IID list for this strata
    iid_list = set(strata_samples['IID'].unique())
    
    # Filter sample list to only include samples in this strata
    filtered_samples = samples_df[samples_df['IID'].isin(iid_list)]
    
    # Save filtered sample list (named to match EXPORT_PLINK's expected pattern)
    filtered_samples.to_csv("${strata}_filtered.pca.harmonized.tsv", sep="\\t", index=False)
    
    print(f"Filtered ${strata}: {len(filtered_samples)} samples retained from {len(samples_df)} total")
    """
}

// ==================================================================================
// MAIN WORKFLOW
// ==================================================================================
workflow {
    
    // ================== Validate Inputs ==================
    if (!params.focus_plink_input || !params.focus_covar_file || !params.focus_pheno_name) {
        error "Missing required parameters: focus_plink_input, focus_covar_file, focus_pheno_name"
    }
    
    // ================== Setup Directories ==================
    def OUTPUT_DIR = "${params.project_dir}/analyses/${params.analysis_name}"
    
    // ================== Load Focus Plink Files ==================
    def focus_prefix = params.focus_plink_input.replaceFirst(/\.pgen$/, '')
    def plen = file(focus_prefix + ".pgen")
    def pvar = file(focus_prefix + ".pvar")
    def psam = file(focus_prefix + ".psam")
    
    plink_files = Channel.of([plen, pvar, psam]).collect()
    
    // ================== Load Covariate File ==================
    covar_file_raw = Channel.fromPath(params.focus_covar_file, checkIfExists: true)
    
    // Rename to match expected pattern for HARMONIZE_CATEGORICAL_COVARS
    RENAME_COVAR_FILE(covar_file_raw)
    covar_file = RENAME_COVAR_FILE.out
    
    // ================== Step 1: Harmonize Categorical Covariates when some cats are too few ====
    HARMONIZE_CATEGORICAL_COVARS(covar_file)
    harmonized_covar = HARMONIZE_CATEGORICAL_COVARS.out
    
    // ================== Step 2: Stratification setup ==================
    if (params.focus_strata_file) {
        strata_file = Channel.fromPath(params.focus_strata_file, checkIfExists: true)
        strata_channel = strata_file.map { f ->
            f.readLines().drop(1).collect { it.split('\t')[2] }.unique().sort()
        }.flatten()

        if (params.longitudinal_flag || params.survival_flag) {
            // For longitudinal/survival: filter harmonized covar by strata now,
            // because GWASGALLOP/GWASCPH need a per-strata samplelist
            FILTER_BY_STRATA(harmonized_covar, strata_file, strata_channel)
            sample_lists = FILTER_BY_STRATA.out
        }
        // For GLM: filtering happens in Step 4 on the pheno.tsv (which has phenotype column)
    } else {
        // No stratification: single analysis group
        sample_lists = harmonized_covar.map { f ->
            def name = f.baseName.replaceAll(/_filtered\.pca\.harmonized\.tsv$/, '')
            tuple(name, f)
        }
    }
    
    // ================== Step 3: Data Export (happens once) ==================
    if (params.longitudinal_flag || params.survival_flag) {
        // For longitudinal/survival: normalize variant IDs to chr:pos:ref:alt first,
        // because survival.R and gallop.py parse IDs by splitting on ':'
        def fileTag = file(focus_prefix).name
        NORMALIZE_VARIANT_IDS(Channel.of(tuple(fileTag, plen, pvar, psam)))

        // RAWFILE_EXPORT expects tuple(fileTag, [pgen, pvar, psam]) — files as a list
        // Use .first() to create a value channel so it broadcasts to all strata samplelists
        norm_plink = NORMALIZE_VARIANT_IDS.out
            .map { normTag, normPgen, normPvar, normPsam -> tuple(normTag, [normPgen, normPvar, normPsam]) }
            .first()

        // Use strata-filtered sample_lists (if stratified), otherwise harmonized_covar
        def export_samplelist
        if (params.focus_strata_file) {
            export_samplelist = sample_lists.map { strata, sl -> sl }
        } else {
            export_samplelist = harmonized_covar
        }
        RAWFILE_EXPORT(norm_plink, export_samplelist)
        export_data = RAWFILE_EXPORT.out.gwas_rawfile.transpose()
        
    } else {
        // For GLM: EXPORT_PLINK runs ONCE on the full harmonized covar.
        // Per-strata filtering is applied to its pheno.tsv output in step 4,
        // so every GWASGLM job gets a file with both phenotype AND covariate columns.
        EXPORT_PLINK(harmonized_covar, params.phenofile)
        export_pheno       = EXPORT_PLINK.out[0]           // A_filtered.pca.pheno.tsv (all samples)
        export_covar_names = EXPORT_PLINK.out[1].first()   // broadcast as value channel
        export_n_covar     = EXPORT_PLINK.out[2].first()   // broadcast as value channel
    }
    
    // ================== Step 4: GWAS Analysis by Strata ==================
    if (params.longitudinal_flag) {
        // Longitudinal analysis: GWASGALLOP — matches main.nf pattern
        GWASGALLOP(export_data, params.phenofile, Channel.value(params.focus_pheno_name))
        gwas_results = GWASGALLOP.out
        
    } else if (params.survival_flag) {
        // Survival analysis: GWASCPH — matches main.nf pattern
        // Pass params.phenofile directly (not as queue channel) so it broadcasts to all chunks
        GWASCPH(export_data, params.phenofile, Channel.value(params.focus_pheno_name))
        gwas_results = GWASCPH.out
        
    } else {
        // Cross-sectional: GWASGLM — one job per strata
        def fileTag = file(focus_prefix).name

        // Filter the pheno.tsv (has phenotype + covariates) per strata — one GWASGLM job each
        if (params.focus_strata_file) {
            filtered_phenos = FILTER_BY_STRATA(export_pheno, strata_file, strata_channel)
        } else {
            filtered_phenos = export_pheno.map { f ->
                def name = f.name.replaceAll(/_filtered\.pca\.pheno\.tsv$/, '')
                tuple(name, f)
            }
        }

        gwas_input = filtered_phenos
            .combine(export_covar_names)
            .combine(export_n_covar)
            .map { strata, pheno, covar_names_f, n_covar_f ->
                tuple(fileTag, plen, plen, psam, pvar, strata, pheno, covar_names_f, n_covar_f)
            }

        GWASGLM(gwas_input, Channel.value(params.focus_pheno_name))
        gwas_results = GWASGLM.out[0].flatten()
    }
    
    // ================== Publish Results ==================
    gwas_results.subscribe { result ->
        if (result instanceof java.nio.file.Path) {
            println "Result: ${result.name}"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
