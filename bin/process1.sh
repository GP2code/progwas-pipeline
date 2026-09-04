#!/bin/bash

# Take the vcf file (genotyped or imputed)
# Resource folder needed and files will be downloaded on-demand if not present
# Process:
# 1. Filter PASS (& R2)
# 2. Liftover to hg38 (if needed)
# 3. Split multiallelics
# 4. Left-normalize on hg38
# 5. Keep only well-behaved SNPs (MAC≥2, ref-aligned, no provisional, no dups, geno<0.1)
# 6. Standard naming chr:pos:ref:alt

## e.g. process1.sh 2 '/data/CARD/PD/imputed_data/CORIELL/chr21.dose.vcf.gz' 0.3 hg19 chr21_cor
## e.g. using docker:
# rm -rf test_output && mkdir -p References test_output && docker run --rm \
#   -v "$PWD:/workspace" \
#   -v "$PWD/References:/workspace/References" \
#   -e RESOURCE_DIR=/workspace/References \
#   -w /workspace \
#   longgwas:slim bash -c "cd test_output && bash ../bin/process1.sh 2 ../example/genotype/chr21.vcf -9 hg19 chr21_test 2>&1 | tee process1.log"


# Parameters
N=$1 # Threads to use
VFILE=$2
R2THRES=$3 # If imputed, give a number for R2 threshold (usually 0.3 - 0.8) -9 otherwise
ASSEMBLY=$4 # [hg18, hg19, hg38]. Define if the liftover is required or not
FILE=$5 # Base file name. can be anything as long as unique

# Resources (References always mounted from host)
RESOURCE_DIR=${RESOURCE_DIR:-./References}
FA=${RESOURCE_DIR}/Genome/hg38.fa.gz
LIFTOVERCHAIN=${RESOURCE_DIR}/liftOver/${ASSEMBLY}ToHg38.over.chain.gz

# Verify reference genomes are available (should be pre-downloaded by DOWNLOAD_REFERENCES process)
if [ ! -f "$FA" ] || [ ! -f "${FA}.fai" ]; then
    echo "ERROR: Reference genome hg38 not found at: $FA"
    echo "Expected files:"
    echo "  - $FA"
    echo "  - ${FA}.fai"
    echo ""
    echo "References should be downloaded by the DOWNLOAD_REFERENCES process before this step."
    echo "If running manually, execute: bin/download_references.sh $ASSEMBLY"
    exit 1
fi

if [ "$ASSEMBLY" != "hg38" ]; then
    ASSEMBLY_FA=${RESOURCE_DIR}/Genome/${ASSEMBLY}.fa.gz
    if [ ! -f "$ASSEMBLY_FA" ] || [ ! -f "${ASSEMBLY_FA}.fai" ] || [ ! -f "$LIFTOVERCHAIN" ]; then
        echo "ERROR: Source assembly ${ASSEMBLY} references not found"
        echo "Expected files:"
        echo "  - $ASSEMBLY_FA"
        echo "  - ${ASSEMBLY_FA}.fai"
        echo "  - $LIFTOVERCHAIN"
        echo ""
        echo "References should be downloaded by the DOWNLOAD_REFERENCES process before this step."
        echo "If running manually, execute: bin/download_references.sh $ASSEMBLY"
        exit 1
    fi
fi

######## start processing ###############################
# Step 1: Filter PASS (&R2 if imputed) and add "chr" prefix if missing
## Different pipeline for imputed and genotyped (R2THRES>0....imputed)
if (( $(awk -v r2="$R2THRES" 'BEGIN {print (r2 > 0)}') ))
then 
    bcftools view -f '.,PASS' \
                  -i "INFO/R2>${R2THRES}" ${VFILE} \
                  -Oz -o ${FILE}_filtered_temp.vcf.gz --threads ${N} # get PASS(or .) variant and R2>R2THRES
else
    bcftools view -f '.,PASS' ${VFILE} \
                  -Oz -o ${FILE}_filtered_temp.vcf.gz --threads ${N}
fi

# Check if chromosomes need "chr" prefix and add if missing
FIRST_CHR=$(bcftools view -H ${FILE}_filtered_temp.vcf.gz | head -1 | cut -f1)
if [[ -z "$FIRST_CHR" ]]; then
    echo "ERROR: No variants found in input VCF after filtering"
    exit 1
fi
if [[ ! "$FIRST_CHR" =~ ^chr ]]; then
    echo "Adding 'chr' prefix to chromosome names..."
    EXPECTED_CHR="chr${FIRST_CHR}"
    bcftools annotate --rename-chrs <(bcftools view -h ${FILE}_filtered_temp.vcf.gz | \
        grep "^##contig" | sed 's/.*ID=\([^,]*\).*/\1/' | \
        awk '{if ($1 !~ /^chr/ && $1 ~ /^[0-9XYM]/) print $1"\tchr"$1; else print $1"\t"$1}') \
        ${FILE}_filtered_temp.vcf.gz -Oz -o ${FILE}_filtered.vcf.gz --threads ${N}
    rm ${FILE}_filtered_temp.vcf.gz
else
    echo "Chromosome names already have 'chr' prefix"
    EXPECTED_CHR="${FIRST_CHR}"
    mv ${FILE}_filtered_temp.vcf.gz ${FILE}_filtered.vcf.gz
fi

# Step 2: Liftover to hg38 if needed (IN VCF FORMAT, BEFORE splitting)
if [[ $ASSEMBLY = hg38 ]]
then
    # Already hg38, just copy
    cp ${FILE}_filtered.vcf.gz ${FILE}_hg38.vcf.gz
else
    # Liftover using bcftools +liftover plugin
    # Note: bcftools +liftover requires both source and destination reference genomes
    # ASSEMBLY_FA already set in reference checking section above
    
    bcftools +liftover ${FILE}_filtered.vcf.gz \
             --threads ${N} \
             -Oz -o ${FILE}_hg38.vcf.gz \
             -- \
             -s ${ASSEMBLY_FA} \
             -f ${FA} \
             -c ${LIFTOVERCHAIN} \
             --reject ${FILE}_rejected.vcf.gz \
             --reject-type z
fi

# Step 3: Split multiallelic variants
bcftools norm -m-both ${FILE}_hg38.vcf.gz \
          -Oz -o ${FILE}_split.vcf.gz --threads ${N}

# Step 4: Convert to plink format (after liftover and split)
if (( $(awk -v r2="$R2THRES" 'BEGIN {print (r2 > 0)}') ))
then
    plink2 --threads ${N} \
           --vcf ${FILE}_split.vcf.gz dosage=DS \
           --make-pgen --allow-extra-chr --autosome-par --out ${FILE}_split
    EXIT_CODE=$?
    if [ $EXIT_CODE -eq 137 ] || [ $EXIT_CODE -eq 143 ]; then
        echo "ERROR: plink2 was killed (likely OOM). Exit code: $EXIT_CODE" >&2
        echo "OOM_ERROR" > ${FILE}_error_type.txt
        exit $EXIT_CODE
    elif [ $EXIT_CODE -ne 0 ]; then
        echo "ERROR: plink2 failed with exit code: $EXIT_CODE" >&2
        exit $EXIT_CODE
    fi
else
    plink2 --threads ${N} \
           --vcf ${FILE}_split.vcf.gz --make-pgen --allow-extra-chr --autosome-par --out ${FILE}_split
    EXIT_CODE=$?
    if [ $EXIT_CODE -eq 137 ] || [ $EXIT_CODE -eq 143 ]; then
        echo "ERROR: plink2 was killed (likely OOM). Exit code: $EXIT_CODE" >&2
        echo "OOM_ERROR" > ${FILE}_error_type.txt
        exit $EXIT_CODE
    elif [ $EXIT_CODE -ne 0 ]; then
        echo "ERROR: plink2 failed with exit code: $EXIT_CODE" >&2
        exit $EXIT_CODE
    fi
fi

# Step 5: Left-normalize using fasta file (hg38)
plink2 --threads ${N} \
       --pfile ${FILE}_split --make-pgen --fa $FA --normalize --sort-vars --out ${FILE}_split_hg38_normalized
EXIT_CODE=$?
if [ $EXIT_CODE -eq 137 ] || [ $EXIT_CODE -eq 143 ]; then
    echo "ERROR: plink2 normalize was killed (likely OOM). Exit code: $EXIT_CODE" >&2
    echo "OOM_ERROR" > ${FILE}_error_type.txt
    exit $EXIT_CODE
elif [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: plink2 normalize failed with exit code: $EXIT_CODE" >&2
    exit $EXIT_CODE
fi

# select relevant snps with more than sigleton (For small cohorts, this process reduces a lot of variants)
# Also filter to keep only the expected chromosome
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized --make-pgen --snps-only just-acgt --mac 2 --chr ${EXPECTED_CHR} --out ${FILE}_split_hg38_normalized_snps
# align ref alt (This will return provisional variants sometime)
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps --make-pgen --ref-from-fa force --fa $FA --out ${FILE}_split_hg38_normalized_snps_temp_aligned
# remove the "Provisional variants" because provisional variants prevents loading on plink2
grep 'PR$' ${FILE}_split_hg38_normalized_snps_temp_aligned.pvar | cut -f3 > ${FILE}_split_hg38_normalized_snps_temp_aligned_provisional.txt || true
# exclude provisional variants first (skip if no provisional variants found)
if [[ -s ${FILE}_split_hg38_normalized_snps_temp_aligned_provisional.txt ]]; then
    plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps --make-pgen --exclude ${FILE}_split_hg38_normalized_snps_temp_aligned_provisional.txt --out ${FILE}_split_hg38_normalized_snps_noprov
else
    # No provisional variants, just copy
    cp ${FILE}_split_hg38_normalized_snps.pgen ${FILE}_split_hg38_normalized_snps_noprov.pgen
    cp ${FILE}_split_hg38_normalized_snps.pvar ${FILE}_split_hg38_normalized_snps_noprov.pvar
    cp ${FILE}_split_hg38_normalized_snps.psam ${FILE}_split_hg38_normalized_snps_noprov.psam
fi
# then align ref/alt independently
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps_noprov --make-pgen --ref-from-fa force --fa $FA --out ${FILE}_split_hg38_normalized_snps_aligned
# remame ID for standard chr:pos:ref:alt
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps_aligned --make-pgen --set-all-var-ids 'chr@:#:$r:$a' --out ${FILE}_split_hg38_normalized_snps_aligned_renamed
# remove dup
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps_aligned_renamed --make-pgen --rm-dup exclude-all --out ${FILE}_split_hg38_normalized_snps_aligned_renamed_uniq
# geno 0.1 filter (lenient because potentially mixed ancestry samples)
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps_aligned_renamed_uniq --make-pgen --geno 0.1 dosage --out ${FILE}_split_hg38_normalized_snps_aligned_renamed_uniq_geno01
# convert to plink pgen format with standardized pvar columns for merging
plink2 --threads ${N} --pfile ${FILE}_split_hg38_normalized_snps_aligned_renamed_uniq_geno01 --make-pgen 'pvar-cols=' --out ${FILE}_p1out