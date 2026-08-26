#!/bin/bash

# Simplified QC pipeline for ancestry-specific data (skip population splitting mode)
# Performs: callrate filter, heterozygosity filter, kinship filter, PCA
# Usage: simple_qc.sh <input_prefix> <kinship_cutoff> <output_prefix> [threads] [mind]

GENO=$1
KINSHIP_CUTOFF=$2
OUT=$3
THREADS=${4:-4}  # Default to 4 threads if not specified
MIND=${5:-0.05}  # Sample call rate threshold (default 0.05 if not provided)

echo "=== Simplified QC Pipeline (skip population splitting) ==="
echo "Input: ${GENO}"
echo "Kinship cutoff: ${KINSHIP_CUTOFF}"
echo "Output: ${OUT}"
echo "Threads: ${THREADS}"
echo "Mind threshold: ${MIND}"

# Step 1: Callrate filtering
echo "Step 1: Callrate filtering..."
plink2 --pfile ${GENO} \
       --mind ${MIND} dosage \
       --make-pgen \
       --threads ${THREADS} \
       --out ${GENO}_callrate

# Step 2: Heterozygosity filtering
echo "Step 2: Heterozygosity filtering..."
plink2 --pfile ${GENO}_callrate \
       --het \
       --threads ${THREADS} \
       --out ${GENO}_callrate_het

# Calculate mean and SD of F coefficient, remove outliers (>3 SD from mean)
awk 'NR>1 {print $1, $5}' ${GENO}_callrate_het.het > ${GENO}_het_values.txt

MEAN=$(awk '{sum+=$2; count++} END {print sum/count}' ${GENO}_het_values.txt)
SD=$(awk -v mean=$MEAN '{sum+=($2-mean)^2; count++} END {print sqrt(sum/count)}' ${GENO}_het_values.txt)
LOWER=$(awk -v mean="$MEAN" -v sd="$SD" 'BEGIN {print mean - 3*sd}')
UPPER=$(awk -v mean="$MEAN" -v sd="$SD" 'BEGIN {print mean + 3*sd}')

echo "Heterozygosity: mean=$MEAN, SD=$SD, range=[$LOWER, $UPPER]"

awk -v lower=$LOWER -v upper=$UPPER 'NR==1 {print "#IID"} NR>1 {if ($5 >= lower && $5 <= upper) print $1}' \
    ${GENO}_callrate_het.het > ${GENO}_het_keep.txt

plink2 --pfile ${GENO}_callrate \
       --keep ${GENO}_het_keep.txt \
       --make-pgen \
       --threads ${THREADS} \
       --out ${GENO}_callrate_het

# Step 3: Kinship filtering
echo "Step 3: Kinship filtering (cutoff=${KINSHIP_CUTOFF})..."
plink2 --pfile ${GENO}_callrate_het \
       --king-cutoff ${KINSHIP_CUTOFF} \
       --make-pgen \
       --threads ${THREADS} \
       --out ${GENO}_callrate_het_king

# Step 4: PCA calculation
echo "Step 4: Computing PCA..."
plink2 --pfile ${GENO}_callrate_het_king \
       --pca approx\
       --threads ${THREADS} \
       --out ${OUT}_pca

# Step 5: Generate kinship table (filter to reduce file size)
# Use lenient threshold (0.0442 = 4th degree relatives) to capture potential relatives
# while reducing file size from GB to MB. Actual filtering done in simple_qc_helper.py
echo "Step 5: Generating kinship table (pairs >= 0.0442, 4th degree relatives)..."
plink2 --pfile ${GENO}_callrate_het_king \
       --make-king-table \
       --king-table-filter 0.0442 \
       --threads ${THREADS} \
       --out ${OUT}_king

# Step 6: Collect outlier information
echo "Step 6: Collecting outlier information..."
TOTAL_SAMPLES=$(wc -l < ${GENO}.psam | awk '{print $1-1}')
AFTER_CALLRATE=$(wc -l < ${GENO}_callrate.psam | awk '{print $1-1}')
AFTER_HET=$(wc -l < ${GENO}_callrate_het.psam | awk '{print $1-1}')
AFTER_KING=$(wc -l < ${GENO}_callrate_het_king.psam | awk '{print $1-1}')

REMOVED_CALLRATE=$((TOTAL_SAMPLES - AFTER_CALLRATE))
REMOVED_HET=$((AFTER_CALLRATE - AFTER_HET))
REMOVED_KING=$((AFTER_HET - AFTER_KING))
TOTAL_REMOVED=$((TOTAL_SAMPLES - AFTER_KING))

cat > ${OUT}_qc_summary.txt << EOF
=== QC Summary ===
Total input samples: ${TOTAL_SAMPLES}
After callrate filter: ${AFTER_CALLRATE} (removed: ${REMOVED_CALLRATE})
After heterozygosity filter: ${AFTER_HET} (removed: ${REMOVED_HET})
After kinship filter: ${AFTER_KING} (removed: ${REMOVED_KING})
Total samples passing QC: ${AFTER_KING}
Total samples removed: ${TOTAL_REMOVED}
EOF

cat ${OUT}_qc_summary.txt

# Step 7: Create HDF5 file and PCA plots
echo "Step 7: Creating HDF5 output and PCA plots..."
python3 $(dirname "$0")/simple_qc_helper.py "${GENO}" "${OUT}"

echo "=== QC Pipeline Complete ==="
