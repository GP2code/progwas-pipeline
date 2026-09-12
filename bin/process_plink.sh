#!/usr/bin/env bash
# Standardize plink files: dedup, optional liftover, normalize to hg38
# Alternative to process1.sh for when input is already in plink format
# Usage: process_plink.sh <threads> <bfile_or_pfile> <assembly> <output_prefix> [r2thres]

set -euo pipefail

N=$1          # Threads
INFILE=$2     # Input bfile/pfile prefix (without .bed/.pgen extension)
ASSEMBLY=$3   # hg19, hg38
OUTPREFIX=$4  # Output file prefix
R2THRES=${5:-0.3}  # R2 threshold (default 0.3 if not provided)

# Resources (References always mounted from host)
RESOURCE_DIR=${RESOURCE_DIR:-/workspace/References}
FA_HG38="${RESOURCE_DIR}/Genome/hg38.fa.gz"
CHAIN_TO38="${RESOURCE_DIR}/liftOver/${ASSEMBLY}ToHg38.over.chain.gz"

# Detect input format (bed or pgen)
if [[ -f "${INFILE}.bed" ]]; then
    INFMT="--bfile"
    echo "[INFO] Input format: PLINK binary (.bed/.bim/.fam)"
elif [[ -f "${INFILE}.pgen" ]]; then
    INFMT="--pfile"
    echo "[INFO] Input format: PLINK2 (.pgen/.pvar/.psam)"
else
    echo "ERROR: No .bed or .pgen found for: ${INFILE}" >&2
    exit 1
fi

# Check reference files
[[ -f "$FA_HG38" ]] || { echo "ERROR: hg38 FASTA not found: $FA_HG38" >&2; exit 1; }
[[ -f "$FA_HG38.fai" ]] || { echo "ERROR: hg38 FASTA index missing: $FA_HG38.fai" >&2; exit 1; }

echo "[INFO] Processing: ${INFILE} -> ${OUTPREFIX}"
echo "[INFO] Assembly: ${ASSEMBLY} -> hg38"
echo "[INFO] Threads: ${N}"

# Step 1: Remove duplicates, keep biallelic SNPs only, MAF>=2
plink2 ${INFMT} "${INFILE}" \
       --set-all-var-ids '@:#:$r:$a' \
       --new-id-max-allele-len 100 truncate \
       --rm-dup exclude-all \
       --snps-only just-acgt \
       --max-alleles 2 \
       --mac 2 \
       --make-pgen \
       --threads "$N" \
       --out "${OUTPREFIX}_dedup"

# Step 1.5: Filter by R2 if INFO field contains R2 scores
if grep -q "R2=" "${OUTPREFIX}_dedup.pvar" 2>/dev/null; then
    echo "[INFO] Filtering variants by R2 >= ${R2THRES}"
    
    NUM_BEFORE=$(grep -vc "^#" "${OUTPREFIX}_dedup.pvar" || echo 0)
    
    # Use PLINK2's native INFO filtering
    plink2 --pfile "${OUTPREFIX}_dedup" \
           --extract-if-info "R2 >= ${R2THRES}" \
           --make-pgen \
           --threads "$N" \
           --out "${OUTPREFIX}_dedup_r2filtered"
    rm -f "${OUTPREFIX}_dedup".{pgen,pvar,psam,log}
    
    NUM_AFTER=$(grep -vc "^#" "${OUTPREFIX}_dedup_r2filtered.pvar" || echo 0)
    echo "[INFO] R2 filter: ${NUM_BEFORE} -> ${NUM_AFTER} variants (removed: $((NUM_BEFORE - NUM_AFTER)))"
    
    WORKPFX="${OUTPREFIX}_dedup_r2filtered"
else
    echo "[INFO] No R2 information found in pvar file, skipping R2 filter"
    WORKPFX="${OUTPREFIX}_dedup"
fi

# Step 2: Liftover if needed
if [[ "$ASSEMBLY" != "hg38" ]]; then
    echo "[INFO] Lifting from ${ASSEMBLY} to hg38"
    [[ -f "$CHAIN_TO38" ]] || { echo "ERROR: Chain file not found: $CHAIN_TO38" >&2; exit 1; }
    
    # Set variant IDs for tracking
    plink2 --pfile "$WORKPFX" \
           --set-all-var-ids '@:#:$r:$a' \
           --new-id-max-allele-len 999 truncate \
           --make-pgen \
           --threads "$N" \
           --out "${OUTPREFIX}_named"
    rm -f "${WORKPFX}".{pgen,pvar,psam,log}
    
    # Prepare BED for liftOver (0-based, chr prefix for UCSC)
    awk 'BEGIN{OFS="\t"} !/^#/{
        chr=$1; pos=$2; id=$3; ref=$4;
        if (chr !~ /^chr/) chr="chr" chr;
        start=pos-1; end=start+length(ref);
        print chr, start, end, id
    }' "${OUTPREFIX}_named.pvar" > "${OUTPREFIX}_lift.in.bed"
    
    # Run liftOver
    liftOver "${OUTPREFIX}_lift.in.bed" "$CHAIN_TO38" \
             "${OUTPREFIX}_lift.out.bed" "${OUTPREFIX}_lift.unmapped"
    
    # Keep first mapping per variant (dedup multi-maps)
    awk '!seen[$4]++' "${OUTPREFIX}_lift.out.bed" > "${OUTPREFIX}_lift.unique.bed"
    
    # Extract kept IDs and new positions
    cut -f4 "${OUTPREFIX}_lift.unique.bed" > "${OUTPREFIX}_keep.ids"
    awk 'BEGIN{OFS="\t"} {print $4, $2+1}' "${OUTPREFIX}_lift.unique.bed" > "${OUTPREFIX}_update.pos"
    awk 'BEGIN{OFS="\t"} {print $4, $1}' "${OUTPREFIX}_lift.unique.bed" > "${OUTPREFIX}_update.chr"
    
    # Apply liftover results
    plink2 --pfile "${OUTPREFIX}_named" \
           --extract "${OUTPREFIX}_keep.ids" \
           --update-map "${OUTPREFIX}_update.pos" \
           --update-chr "${OUTPREFIX}_update.chr" 2 1 \
           --sort-vars \
           --make-pgen \
           --threads "$N" \
           --out "${OUTPREFIX}_hg38"
    rm -f "${OUTPREFIX}_named".{pgen,pvar,psam,log}
    rm -f "${OUTPREFIX}"_lift.* "${OUTPREFIX}"_keep.ids "${OUTPREFIX}"_update.*
    
    WORKPFX="${OUTPREFIX}_hg38"
    
    echo "[INFO] Liftover: $(wc -l < "${OUTPREFIX}_lift.in.bed") -> $(wc -l < "${OUTPREFIX}_lift.unique.bed")"
else
    echo "[INFO] Input already hg38, ensuring chr prefix"
    # Add chr prefix if not present (variant_ID, new_chr)
    awk 'BEGIN{OFS="\t"} /^#/ {next} $1 !~ /^chr/ {print $3, "chr"$1; next} {print $3, $1}' \
        "$WORKPFX.pvar" > "${OUTPREFIX}_chr_update.txt"
    plink2 --pfile "$WORKPFX" \
           --update-chr "${OUTPREFIX}_chr_update.txt" 2 1 \
           --sort-vars \
           --make-pgen \
           --threads "$N" \
           --out "${OUTPREFIX}_hg38"
    rm -f "${WORKPFX}".{pgen,pvar,psam,log} "${OUTPREFIX}_chr_update.txt"
    WORKPFX="${OUTPREFIX}_hg38"
fi

# Step 3: Normalize (left-align, split multiallelic)
echo "[INFO] Normalizing variants (left-align)"
plink2 --pfile "$WORKPFX" \
       --fa "$FA_HG38" \
       --normalize \
       --rm-dup exclude-all \
       --sort-vars \
       --make-pgen \
       --threads "$N" \
       --out "${OUTPREFIX}_norm"
rm -f "${WORKPFX}".{pgen,pvar,psam,log}

# Step 4: geno 0.1 and convert to hard-call
## without this process, raw file has dosage and inconsistent with VCF-based processing)
plink2 --pfile "${OUTPREFIX}_norm" \
       --geno 0.1 \
       --make-pgen \
       --threads "$N" \
       --out "${OUTPREFIX}_geno_hc"
rm -f "${OUTPREFIX}_norm".{pgen,pvar,psam,log}

# Step 5: Align REF/ALT to reference FASTA
echo "[INFO] Aligning REF/ALT to hg38 reference"
plink2 --pfile "${OUTPREFIX}_geno_hc" \
       --fa "$FA_HG38" \
       --ref-from-fa force \
       --make-pgen \
       --threads "$N" \
       --out "${OUTPREFIX}_refalign"
rm -f "${OUTPREFIX}_geno_hc".{pgen,pvar,psam,log}

# Step 6: Set final variant IDs and output as pgen format (compatible with MERGER_CHRS)
plink2 --pfile "${OUTPREFIX}_refalign" \
       --set-all-var-ids 'chr@:#:$r:$a' \
       --new-id-max-allele-len 999 truncate \
       --rm-dup exclude-all \
       --chr 1-22,X,Y \
       --sort-vars \
       --make-pgen \
       --threads "$N" \
       --out "${OUTPREFIX}_p1out"
rm -f "${OUTPREFIX}_refalign".{pgen,pvar,psam,log}

echo "[INFO] Complete: ${OUTPREFIX}_p1out.{pgen,pvar,psam}"
