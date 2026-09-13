// ==============================================================================
// GWAS Execution Module
// ==============================================================================
// Consolidated module containing GWAS analysis processes:
// - GWASGLM: Generalized Linear Model analysis (cross-sectional)
// - GWASGALLOP: Longitudinal analysis using GALLOP
// - GWASCPH: Survival analysis using Cox Proportional Hazards
// ==============================================================================

process GWASGLM {
  scratch true
  label 'medium'

  publishDir "${params.project_dir}/analyses/${params.genetic_cache_key}/${params.analysis_name}/results", mode: 'copy', overwrite: true

  input:
    tuple val(fileTag), path(plog), path(pgen), path(psam), path(pvar), val(pop_studyarm), path(samplelist), path(covar_names_file), path(n_covar_file)
    val phenonames

  output:
    path "*.results"
    path "${outfile}.manifest.tsv"
    // manifest maps each result file to its key (pop_studyarm_phenotype).
    // Filename must be unique per task since all chromosomes share one publishDir.

  script:
    // No `def` here -- output: needs outfile in the shared process scope.
    outfile = "${pop_studyarm}_${fileTag}"
    // Convert phenonames to space-separated string for plink2
    // Handle both String and List input formats
    def pheno_list = phenonames instanceof List ? phenonames.join(' ') : phenonames.toString().replaceAll(/[\[\]'"]/, '').trim()
    // Get basename for plink files (without extension)
    def pfile_base = pgen.getSimpleName()

    """
    set -x

    # Read pre-computed covariate names and count from EXPORT_PLINK
    COVAR_NAMES=\$(cat ${covar_names_file})
    N_COVAR=\$(cat ${n_covar_file})
    
    echo "Using covariates: \${COVAR_NAMES}"
    echo "Total covariates: \${N_COVAR}"
    echo "Processing phenotypes: ${pheno_list}"

    awk -F'\t' 'NR==1{n=NF; next} NF!=n{print "ERROR: ${samplelist} line "NR" has "NF" fields, expected "n" (from header). Line (truncated): "substr(\$0,1,200); exit 1}' "${samplelist}"

    # Build FID/IID keep list from phenotype/covariate table using FID lookup from psam.
    make_keep_iid.py --input "${samplelist}" --output "${outfile}.keep.iid.tsv" --psam "${psam}"
    
    # Note: ${samplelist} contains all samples with standardized covariates from EXPORT_PLINK
    # plink2 --glm automatically excludes samples with missing phenotype values per phenotype
    # Passing multiple phenotypes is much more efficient than iterating

    # Build plink2 command with optional interaction parameters
    if [ -n "${params.covar_interact}" ]; then
        # With interaction: test SNP main effect and SNP*covariate interaction.
        #
        # plink2 numbers --glm predictors as:
        #   1                      ADD (the genotype)
        #   2 .. N_COVAR+1         the covariates
        #   N_COVAR+2 .. 2*N_COVAR+1  the ADDxcovariate interaction terms
        #
        # Crucially, the covariates are ordered by their COLUMN ORDER IN THE
        # COVARIATE FILE, not by the order given to --covar-name. So the index of
        # the ADDx<covar_interact> term depends on where that covariate sits in
        # the file, and hardcoding N_COVAR+2 silently selects the interaction
        # with whichever covariate happens to come first.
        INTERACT_RANK=\$(awk -F'\t' -v names="\${COVAR_NAMES}" -v target="${params.covar_interact}" '
            NR==1 {
                n = split(names, want, ",")
                for (j = 1; j <= n; j++) sel[want[j]] = 1
                k = 0
                for (i = 1; i <= NF; i++) {
                    if (\$i in sel) {
                        k++
                        if (\$i == target) { print k; exit 0 }
                    }
                }
                exit 1
            }' "${samplelist}")

        if [ -z "\${INTERACT_RANK}" ]; then
            echo "ERROR: interaction covariate '${params.covar_interact}' is not among the model covariates (\${COVAR_NAMES}) in ${samplelist}" >&2
            exit 1
        fi

        # Index of ADDx<covar_interact> in plink2's ORIGINAL numbering.
        INTERACTION_IDX=\$((N_COVAR + 1 + INTERACT_RANK))
        # --tests indices refer to positions AFTER --parameters filtering, where
        # the retained terms are ADD (1), the N_COVAR covariates (2..N_COVAR+1)
        # and the single interaction term (N_COVAR+2).
        TEST_IDX=\$((N_COVAR + 2))
        echo "Interaction covariate '${params.covar_interact}' is covariate #\${INTERACT_RANK} in file order -> parameter \${INTERACTION_IDX}"

        plink2 --pfile ${pfile_base} \
                --glm interaction omit-ref cols=+beta,+a1freq \
                --pheno "${samplelist}" \
                --pheno-name ${pheno_list} \
                --covar "${samplelist}" \
                --covar-name \${COVAR_NAMES} \
                --keep "${outfile}.keep.iid.tsv" \
                --output-chr chrM \
                --maf ${params.minor_allele_freq} \
                --mac ${params.minor_allele_ct} \
                --hwe ${params.hwe} \
                --geno ${params.geno} \
                --parameters 1-\$((N_COVAR + 1)),\${INTERACTION_IDX} \
                --tests 1,\${TEST_IDX} \
                --threads ${task.cpus} \
                --memory ${task.memory.toMega()} \
                --out ${outfile}_all_vars
        
        # Reorganize results: ADD as base, join interaction and 2DF columns
        # Process all phenotype output files
        INTERACT_TEST="ADDx${params.covar_interact}"
        for glm_file in ${outfile}_all_vars.*.glm.{linear,logistic.hybrid}; do
            if [ -f "\${glm_file}" ]; then
                output_file=\${glm_file/_all_vars/}
                awk -v interact_test="\${INTERACT_TEST}" 'BEGIN{FS="\t"; OFS="\t"} 
                     NR==1 {
                         # Find column indices
                         for(i=1;i<=NF;i++) {
                             if(\$i=="ID") idcol=i;
                             if(\$i=="TEST") testcol=i;
                             if(\$i=="BETA") betacol=i;
                             if(\$i=="SE") secol=i;
                             if(\$i=="P") pcol=i;
                         }
                         # Print base header plus interaction columns
                         print \$0, "BETA_INT", "SE_INT", "P_INT", "CORR_INT", "P_2DF", "INTERACTION", "MODEL";
                         next;
                     }
                     {
                         id = \$idcol;
                         test = \$testcol;
                         
                         if(test == "ADD") {
                             # Store base ADD row
                             add[id] = \$0;
                         } else if(test == interact_test) {
                             # Store interaction columns for specific interaction term
                             interact_beta[id] = \$betacol;
                             interact_se[id] = \$secol;
                             interact_p[id] = \$pcol;
                         } else if(test == "USER_2DF") {
                             # Store 2DF p-value
                             twodf_p[id] = \$pcol;
                         }
                     }
                     END {
                         # Output combined rows
                         for(id in add) {
                             beta_i = (id in interact_beta) ? interact_beta[id] : "NA";
                             se_i = (id in interact_se) ? interact_se[id] : "NA";
                             p_i = (id in interact_p) ? interact_p[id] : "NA";
                             p_2df = (id in twodf_p) ? twodf_p[id] : "NA";
                             print add[id], beta_i, se_i, p_i, "NA", p_2df, "${params.covar_interact}", "GLM";
                         }
                     }' "\${glm_file}" > "\${output_file}"
            fi
        done
        
    else
        # Standard analysis without interaction
        plink2 --pfile ${pfile_base} \
                --glm hide-covar omit-ref cols=+beta,+a1freq \
                --pheno "${samplelist}" \
                --pheno-name ${pheno_list} \
                --covar "${samplelist}" \
                --covar-name \${COVAR_NAMES} \
                --keep "${outfile}.keep.iid.tsv" \
                --output-chr chrM \
                --maf ${params.minor_allele_freq} \
                --mac ${params.minor_allele_ct} \
                --hwe ${params.hwe} \
                --geno ${params.geno} \
                --threads ${task.cpus} \
                --memory ${task.memory.toMega()} \
                --out ${outfile}
    fi

    # Rename all phenotype output files to .results extension and create manifest
    echo -e "key\tfilename" > ${outfile}.manifest.tsv
    for result_file in ${outfile}.*.glm.{logistic.hybrid,linear}; do
        if [ -f "\${result_file}" ]; then
            new_name="\${result_file%.glm.*}.results"
            mv "\${result_file}" "\${new_name}"
            # In the interaction branch the awk above already appended the
            # BETA_INT..MODEL block, so only pad here when it did not run --
            # otherwise those seven columns end up in the file twice.
            if [ -z "${params.covar_interact}" ]; then
                awk 'BEGIN{OFS="\t"} NR==1{print \$0,"BETA_INT","SE_INT","P_INT","CORR_INT","P_2DF","INTERACTION","MODEL"; next} {print \$0,"NA","NA","NA","NA","NA","NA","GLM"}' "\${new_name}" > "\${new_name}.tmp" && mv "\${new_name}.tmp" "\${new_name}"
            fi

            # Extract phenotype name from filename
            # Pattern: pop_studyarm_fileTag.phenotype.results
            phenotype=\$(basename "\${new_name}" .results | rev | cut -d'.' -f1 | rev)
            key="${pop_studyarm}_\${phenotype}"
            echo -e "\${key}\t\${new_name}" >> ${outfile}.manifest.tsv
        fi
    done
    """
}

process GWASGALLOP {
  scratch true
  label 'medium'

  publishDir "${params.project_dir}/analyses/${params.genetic_cache_key}/${params.analysis_name}/results", mode: 'copy', overwrite: true, enabled: params.publish_gwas_results

  input:
    tuple val(fileTag), path(samplelist), path(rawfile)
    path x, stageAs: 'phenotypes.tsv'
    each phenoname

  output:
    tuple env(KEY), path("*.gallop")

  script:
    def m = []
    def outfile = rawfile.getName()
    m = outfile =~ /(.+)\.raw$/
    outfile = "${m[0][1]}"

    def getkey = []
    def pop_studyarm = samplelist.getName()
    getkey = pop_studyarm =~ /(.+)_filtered\.pca(?:\.harmonized)?\.tsv$/
    pop_studyarm = getkey[0][1]

    """
    set -x
    KEY="${pop_studyarm}_${phenoname}"

    gallop.py --gallop \
           --rawfile ${rawfile} \
           --pheno-file "phenotypes.tsv" \
           --pheno-name "${phenoname}" \
           --covar-file ${samplelist} \
           --covar-numeric ${params.covar_numeric} \
           ${params.covar_categorical ? "--covar-categorical ${params.covar_categorical}" : ""} \
           --time-name ${params.time_col} \
           --out "${outfile}"
    """
}

process GWASCPH {
  scratch true
  label 'small'

  publishDir "${params.project_dir}/analyses/${params.genetic_cache_key}/${params.analysis_name}/results", mode: 'copy', overwrite: true, enabled: params.publish_gwas_results

  input:
    tuple val(fileTag), path(samplelist), path(rawfile)
    path x, stageAs: 'phenotypes.tsv'
    each phenoname

  output:
    tuple env(KEY), path("*.coxph")
  
  script:
    def m = []
    def outfile = rawfile.getName()
    m = outfile =~ /(.+)\.raw$/
    outfile = "${m[0][1]}.coxph"

    def getkey = []
    def pop_studyarm = samplelist.getName()
    getkey = pop_studyarm =~ /(.+)_filtered\.pca(?:\.harmonized)?\.tsv$/
    pop_studyarm = getkey[0][1]

    """
    set -x
    KEY="${pop_studyarm}_${phenoname}"
    
    echo "Processing: ${rawfile.name}"
    echo "Available files:"
    ls -la *.raw *.tsv 2>/dev/null || echo "No files found"
    
    survival.R --rawfile ${rawfile} \
               --pheno-file "phenotypes.tsv" \
               --covar-file ${samplelist} \
               --covar-numeric "${params.covar_numeric}" \
               --covar-categorical "${params.covar_categorical}" \
               ${params.covar_interact ? "--covar-interact \"${params.covar_interact}\"" : ""} \
               --pheno-name "${phenoname}" \
               --time-col "${params.time_col}" \
               --out ${outfile}
    """
}
