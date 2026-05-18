#!/bin/bash

#SBATCH --account=rrg-vmooser
#SBATCH --job-name=03_intersect_LDprune
#SBATCH --output=03_intersect_LDprune.out
#SBATCH --error=03_intersect_LDprune.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

# --------------------------------------------------------------
# Script efficiency (61143707)
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 4
# CPU Utilized: 01:20:55
# CPU Efficiency: 24.21% of 05:34:12 core-walltime
# Job Wall-clock time: 01:23:33
# Memory Utilized: 14.44 GB
# Memory Efficiency: 45.13% of 32.00 GB (32.00 GB/node)

# --------------------------------------------------------------
# Load necessary modules
# --------------------------------------------------------------
module purge
module load StdEnv/2023
module load gcc/12.3
module load bcftools/1.19
module load plink/2.0.0-a.6.32

# Set threads for parallel processing
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# --------------------------------------------------------------
# Input directories/files
# --------------------------------------------------------------
# QCed CARTaGENE VCFs from script 01
CAG_QC_VCF_DIR="/lustre07/scratch/chanalex/CARTaGENE/PCA_QC/04_qc_vcf_for_reference_intersection"

# QCed HGDP-1000G VCFs from script 02
REF_QC_VCF_DIR="/lustre07/scratch/chanalex/HGDP-1KG/PCA_QC/05_qc_vcf_for_CARTaGENE_intersection"

# High-LD regions to exclude before LD pruning 
# (wget from University of Michigan on 2026/05/18)
HIGH_LD_RANGES="/lustre06/project/6061810/chanalex/phd_thesis/Data/high-LD-regions-hg38-GRCh38.txt"

# --------------------------------------------------------------
# Output directories
# --------------------------------------------------------------
BASE_OUT="/lustre07/scratch/chanalex/CARTaGENE_HGDP-1KG/PCA_projection"

ISEC_WORK_DIR="${BASE_OUT}/01_isec_work"

REF_SHARED_BCF_DIR="${BASE_OUT}/02_HGDP_1KG_shared_bcf"
CAG_SHARED_BCF_DIR="${BASE_OUT}/03_CARTaGENE_shared_bcf"

REF_SHARED_PGEN_DIR="${BASE_OUT}/04_HGDP_1KG_shared_pgen"
CAG_SHARED_PGEN_DIR="${BASE_OUT}/05_CARTaGENE_shared_pgen"

LD_DIR="${BASE_OUT}/06_reference_LD_pruning"

REF_LDPRUNED_PGEN_DIR="${BASE_OUT}/07_HGDP_1KG_LDpruned_pgen"
CAG_LDPRUNED_PGEN_DIR="${BASE_OUT}/08_CARTaGENE_LDpruned_pgen"

MERGED_DIR="${BASE_OUT}/09_merged_autosomes"

LOG_DIR="${BASE_OUT}/logs"

mkdir -p \
  "$ISEC_WORK_DIR" \
  "$REF_SHARED_BCF_DIR" \
  "$CAG_SHARED_BCF_DIR" \
  "$REF_SHARED_PGEN_DIR" \
  "$CAG_SHARED_PGEN_DIR" \
  "$LD_DIR" \
  "$REF_LDPRUNED_PGEN_DIR" \
  "$CAG_LDPRUNED_PGEN_DIR" \
  "$MERGED_DIR" \
  "$LOG_DIR"

# --------------------------------------------------------------
# Parameters
# --------------------------------------------------------------
CHROMOSOMES=($(seq -f "chr%g" 1 22))

# LD pruning parameters
# Parameter settings based on: https://www.cog-genomics.org/plink/2.0/ld
# LD step defaults to 1 variant (i.e., window moves by 1 variant at a time)
# when using a kilobase-based window
LD_WINDOW="200kb"
LD_R2=0.5

# --------------------------------------------------------------
# Step 1: Exact CHROM:POS:REF:ALT intersection
#
# -n=2    retain variants present in both datasets
# -c none require exact REF/ALT match
# -w1     write records from first input: HGDP-1000G
# -w2     write records from second input: CARTaGENE
#
# Output BCF rather than VCF for speed/size.
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do

  REF_VCF="${REF_QC_VCF_DIR}/${CHR}.HGDP_1KG.QC.vcf.gz"
  CAG_VCF="${CAG_QC_VCF_DIR}/${CHR}.CARTaGENE.QC.vcf.gz"

  REF_ISEC_DIR="${ISEC_WORK_DIR}/${CHR}.ref_isec"
  CAG_ISEC_DIR="${ISEC_WORK_DIR}/${CHR}.cag_isec"

  REF_OUT_BCF="${REF_SHARED_BCF_DIR}/${CHR}.HGDP_1KG.QC.shared_with_CARTaGENE.bcf"
  CAG_OUT_BCF="${CAG_SHARED_BCF_DIR}/${CHR}.CARTaGENE.QC.shared_with_HGDP_1KG.bcf"

  LOG_FILE="${LOG_DIR}/${CHR}.step1.exact_intersection.log"

  echo "[$(date)] Step 1: Exact HGDP-1000G-CARTaGENE intersection for $CHR"

  if [[ ! -f "$REF_VCF" ]]; then
    echo "ERROR: Missing HGDP-1000G QC VCF: $REF_VCF" >&2
    exit 1
  fi

  if [[ ! -f "$CAG_VCF" ]]; then
    echo "ERROR: Missing CARTaGENE QC VCF: $CAG_VCF" >&2
    exit 1
  fi

  if [[ ! -f "${REF_VCF}.tbi" && ! -f "${REF_VCF}.csi" ]]; then
    bcftools index -t "$REF_VCF"
  fi

  if [[ ! -f "${CAG_VCF}.tbi" && ! -f "${CAG_VCF}.csi" ]]; then
    bcftools index -t "$CAG_VCF"
  fi

  rm -rf "$REF_ISEC_DIR" "$CAG_ISEC_DIR"
  mkdir -p "$REF_ISEC_DIR" "$CAG_ISEC_DIR"

  # Write HGDP-1000G-side overlapping records
  bcftools isec \
    --threads "$THREADS" \
    -n=2 \
    -c none \
    -w1 \
    -Ob \
    -p "$REF_ISEC_DIR" \
    "$REF_VCF" \
    "$CAG_VCF" \
    > "$LOG_FILE" 2>&1

  rm -f "$REF_OUT_BCF" "$REF_OUT_BCF.csi"
  mv "${REF_ISEC_DIR}/0000.bcf" "$REF_OUT_BCF"
  bcftools index "$REF_OUT_BCF"

  # Write CARTaGENE-side overlapping records
  bcftools isec \
    --threads "$THREADS" \
    -n=2 \
    -c none \
    -w2 \
    -Ob \
    -p "$CAG_ISEC_DIR" \
    "$REF_VCF" \
    "$CAG_VCF" \
    >> "$LOG_FILE" 2>&1

  rm -f "$CAG_OUT_BCF" "$CAG_OUT_BCF.csi"
  mv "${CAG_ISEC_DIR}/0001.bcf" "$CAG_OUT_BCF"
  bcftools index "$CAG_OUT_BCF"

  REF_N=$(bcftools index -n "$REF_OUT_BCF")
  CAG_N=$(bcftools index -n "$CAG_OUT_BCF")

  echo "$CHR shared HGDP-1000G variants: $REF_N" | tee -a "$LOG_FILE"
  echo "$CHR shared CARTaGENE variants: $CAG_N" | tee -a "$LOG_FILE"

  if [[ "$REF_N" -eq 0 || "$CAG_N" -eq 0 ]]; then
    echo "ERROR: Zero shared variants for $CHR. Check chromosome naming and REF/ALT matching." >&2
    exit 1
  fi

done

# --------------------------------------------------------------
# Step 2: Convert shared BCFs to PLINK2 PGEN
#         Use the same CHR:POS:REF:ALT variant IDs in both datasets.
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do

  REF_IN_BCF="${REF_SHARED_BCF_DIR}/${CHR}.HGDP_1KG.QC.shared_with_CARTaGENE.bcf"
  CAG_IN_BCF="${CAG_SHARED_BCF_DIR}/${CHR}.CARTaGENE.QC.shared_with_HGDP_1KG.bcf"

  REF_OUT_PREFIX="${REF_SHARED_PGEN_DIR}/${CHR}.HGDP_1KG.QC.shared_with_CARTaGENE"
  CAG_OUT_PREFIX="${CAG_SHARED_PGEN_DIR}/${CHR}.CARTaGENE.QC.shared_with_HGDP_1KG"

  LOG_FILE="${LOG_DIR}/${CHR}.step2.convert_shared_to_pgen.log"

  echo "[$(date)] Step 2: Converting shared variants to PGEN for $CHR"

  plink2 \
    --bcf "$REF_IN_BCF" \
    --const-fid 0 \
    --snps-only just-acgt \
    --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' \
    --rm-dup exclude-all \
    --make-pgen \
    --threads "$THREADS" \
    --out "$REF_OUT_PREFIX" \
    > "$LOG_FILE" 2>&1

  plink2 \
    --bcf "$CAG_IN_BCF" \
    --const-fid 0 \
    --snps-only just-acgt \
    --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' \
    --rm-dup exclude-all \
    --make-pgen \
    --threads "$THREADS" \
    --out "$CAG_OUT_PREFIX" \
    >> "$LOG_FILE" 2>&1

done

# --------------------------------------------------------------
# Step 3: LD prune shared variants in HGDP-1000G reference
#
# Since HGDP-1000G will define the PCA axes, LD pruning is done
# in the reference panel. The resulting SNP list is then applied
# to CARTaGENE.
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do

  REF_IN_PREFIX="${REF_SHARED_PGEN_DIR}/${CHR}.HGDP_1KG.QC.shared_with_CARTaGENE"
  LD_OUT_PREFIX="${LD_DIR}/${CHR}.HGDP_1KG.reference_LD"

  LOG_FILE="${LOG_DIR}/${CHR}.step3.reference_LD_prune.log"

  echo "[$(date)] Step 3: LD pruning in HGDP-1000G for $CHR"

  plink2 \
    --pfile "$REF_IN_PREFIX" \
    --exclude range "$HIGH_LD_RANGES" \
    --indep-pairwise "$LD_WINDOW" "$LD_R2" \
    --threads "$THREADS" \
    --out "$LD_OUT_PREFIX" \
    > "$LOG_FILE" 2>&1

done

cat "${LD_DIR}"/chr*.HGDP_1KG.reference_LD.prune.in \
  | sort -u \
  > "${LD_DIR}/HGDP_1KG.reference_LDpruned_variant_ids.txt"

echo "Total reference LD-pruned variants:"
wc -l "${LD_DIR}/HGDP_1KG.reference_LDpruned_variant_ids.txt"

# --------------------------------------------------------------
# Step 4: Apply the reference-derived LD-pruned SNP list to both
#         HGDP-1000G and CARTaGENE.
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do

  REF_IN_PREFIX="${REF_SHARED_PGEN_DIR}/${CHR}.HGDP_1KG.QC.shared_with_CARTaGENE"
  CAG_IN_PREFIX="${CAG_SHARED_PGEN_DIR}/${CHR}.CARTaGENE.QC.shared_with_HGDP_1KG"

  REF_OUT_PREFIX="${REF_LDPRUNED_PGEN_DIR}/${CHR}.HGDP_1KG.QC.shared_LDpruned"
  CAG_OUT_PREFIX="${CAG_LDPRUNED_PGEN_DIR}/${CHR}.CARTaGENE.QC.shared_LDpruned"

  LOG_FILE="${LOG_DIR}/${CHR}.step4.apply_LDpruned_list.log"

  echo "[$(date)] Step 4: Applying LD-pruned SNP list for $CHR"

  plink2 \
    --pfile "$REF_IN_PREFIX" \
    --extract "${LD_DIR}/${CHR}.HGDP_1KG.reference_LD.prune.in" \
    --make-pgen \
    --threads "$THREADS" \
    --out "$REF_OUT_PREFIX" \
    > "$LOG_FILE" 2>&1

  plink2 \
    --pfile "$CAG_IN_PREFIX" \
    --extract "${LD_DIR}/${CHR}.HGDP_1KG.reference_LD.prune.in" \
    --make-pgen \
    --threads "$THREADS" \
    --out "$CAG_OUT_PREFIX" \
    >> "$LOG_FILE" 2>&1

done

# --------------------------------------------------------------
# Step 5: Merge autosomes into one HGDP-1000G file and one
#         CARTaGENE file.
#
# These merged files will be used in the next script:
# reference PCA + CARTaGENE projection.
# --------------------------------------------------------------
REF_PMERGE_LIST="${MERGED_DIR}/HGDP_1KG_LDpruned_pmerge_list.txt"
CAG_PMERGE_LIST="${MERGED_DIR}/CARTaGENE_LDpruned_pmerge_list.txt"

: > "$REF_PMERGE_LIST"
: > "$CAG_PMERGE_LIST"

for CHR in "${CHROMOSOMES[@]}"; do
  echo "${REF_LDPRUNED_PGEN_DIR}/${CHR}.HGDP_1KG.QC.shared_LDpruned" >> "$REF_PMERGE_LIST"
  echo "${CAG_LDPRUNED_PGEN_DIR}/${CHR}.CARTaGENE.QC.shared_LDpruned" >> "$CAG_PMERGE_LIST"
done

REF_MERGED_PREFIX="${MERGED_DIR}/HGDP_1KG.QC.shared_LDpruned.autosomes"
CAG_MERGED_PREFIX="${MERGED_DIR}/CARTaGENE.QC.shared_LDpruned.autosomes"

echo "[$(date)] Step 5: Merging HGDP-1000G autosomes"

plink2 \
  --pmerge-list "$REF_PMERGE_LIST" \
  --make-pgen \
  --threads "$THREADS" \
  --out "$REF_MERGED_PREFIX" \
  > "${LOG_DIR}/step5.merge_HGDP_1KG_autosomes.log" 2>&1

echo "[$(date)] Step 5: Merging CARTaGENE autosomes"

plink2 \
  --pmerge-list "$CAG_PMERGE_LIST" \
  --make-pgen \
  --threads "$THREADS" \
  --out "$CAG_MERGED_PREFIX" \
  > "${LOG_DIR}/step5.merge_CARTaGENE_autosomes.log" 2>&1

echo "[$(date)] Shared-variant intersection and reference LD pruning complete."
echo "Reference LD-pruned merged PGEN:"
echo "  $REF_MERGED_PREFIX"
echo "CARTaGENE LD-pruned merged PGEN:"
echo "  $CAG_MERGED_PREFIX"
echo "LD-pruned SNP list:"
echo "  ${LD_DIR}/HGDP_1KG.reference_LDpruned_variant_ids.txt"