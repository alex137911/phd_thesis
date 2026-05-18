#!/bin/bash

#SBATCH --account=rrg-vmooser
#SBATCH --job-name=04_referencePCA_HGDP_1KG
#SBATCH --output=04_referencePCA_HGDP_1KG.out
#SBATCH --error=04_referencePCA_HGDP_1KG.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

set -euo pipefail

# --------------------------------------------------------------
# Script efficiency (61169769)
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
module load plink/2.0.0-a.6.32

# Set threads for parallel processing
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# --------------------------------------------------------------
# Input directories/files
# --------------------------------------------------------------
REF_MERGED_PREFIX="/lustre07/scratch/chanalex/CARTaGENE_HGDP-1KG/PCA_projection/09_merged_autosomes/HGDP_1KG.QC.shared_LDpruned.autosomes"

# Check for input files
if [[ ! -f "${REF_MERGED_PREFIX}.pgen" ]]; then
  echo "ERROR: Missing reference merged PGEN: ${REF_MERGED_PREFIX}.pgen" >&2
  exit 1
fi

if [[ ! -f "${REF_MERGED_PREFIX}.pvar" ]]; then
  echo "ERROR: Missing reference merged PVAR: ${REF_MERGED_PREFIX}.pvar" >&2
  exit 1
fi

if [[ ! -f "${REF_MERGED_PREFIX}.psam" ]]; then
  echo "ERROR: Missing reference merged PSAM: ${REF_MERGED_PREFIX}.psam" >&2
  exit 1
fi

# --------------------------------------------------------------
# Output directories
# --------------------------------------------------------------
BASE_OUT="/lustre07/scratch/chanalex/CARTaGENE_HGDP-1KG/PCA_projection"
PCA_DIR="${BASE_OUT}/10_reference_PCA"
LOG_DIR="${BASE_OUT}/logs"

mkdir -p "$PCA_DIR" "$KING_DIR" "$LOG_DIR"

N_PCS=20

# Directory to store KING output (relatedness estimates)
KING_PREFIX="${KING_DIR}/HGDP_1KG.reference.king_cutoff_${KING_CUTOFF}"

# Directory to store PCA output (reference PCA and self-projection)
REF_PCA_PREFIX="${PCA_DIR}/HGDP_1KG.reference_PCA"
REF_SELF_PROJECT_PREFIX="${PCA_DIR}/HGDP_1KG.reference_self_projected"

# --------------------------------------------------------------
# Parameters
# --------------------------------------------------------------
# Following: https://www.kingrelatedness.com/manual.shtml
# Remove one sample from each related pair above this KING kinship threshold.
# 0.177 ≈ remove duplicate/MZ + first-degree relatives.
# 0.0884 ≈ more conservative; removes up to around second-degree relatives.
KING_CUTOFF=0.0884

# Output before and after KING pruning (unrelated samples)
echo "Reference samples before KING pruning:"
awk 'NR > 1 {n++} END {print n}' "${REF_MERGED_PREFIX}.psam"

echo "Reference variants:"
grep -vc '^#' "${REF_MERGED_PREFIX}.pvar"

# --------------------------------------------------------------
# Step 1: KING relatedness pruning in HGDP/1000G reference
#
# This writes:
#   ${KING_PREFIX}.king.cutoff.in.id   = unrelated/reference-training samples
#   ${KING_PREFIX}.king.cutoff.out.id  = excluded related/duplicate samples
# --------------------------------------------------------------
echo "[$(date)] Step 1: Running KING relatedness pruning"

plink2 \
  --pfile "$REF_MERGED_PREFIX" \
  --king-cutoff "$KING_CUTOFF" \
  --threads "$THREADS" \
  --out "$KING_PREFIX" \
  > "${LOG_DIR}/step1.HGDP_1KG.king_cutoff.log" 2>&1

if [[ ! -f "${KING_PREFIX}.king.cutoff.in.id" ]]; then
  echo "ERROR: KING retained-sample file was not created." >&2
  exit 1
fi

echo "Reference samples retained for PCA training:"
grep -vc '^#' "${KING_PREFIX}.king.cutoff.in.id" || wc -l "${KING_PREFIX}.king.cutoff.in.id"

echo "Reference samples excluded by KING pruning:"
if [[ -f "${KING_PREFIX}.king.cutoff.out.id" ]]; then
  grep -vc '^#' "${KING_PREFIX}.king.cutoff.out.id" || wc -l "${KING_PREFIX}.king.cutoff.out.id"
else
  echo "0"
fi

# --------------------------------------------------------------
# Step 2: Compute PCA in unrelated HGDP-1000G reference subset
#
# --freq counts:
#   writes allele counts/frequencies from the unrelated reference subset
#
# --pca allele-wts:
#   writes allele weights for projection
# --------------------------------------------------------------
echo "[$(date)] Step 2: Computing PCA in unrelated HGDP-1000G reference subset"

plink2 \
  --pfile "$REF_MERGED_PREFIX" \
  --keep "${KING_PREFIX}.king.cutoff.in.id" \
  --freq counts \
  --pca "$N_PCS" allele-wts vcols=chrom,ref,alt \
  --threads "$THREADS" \
  --out "$REF_PCA_PREFIX" \
  > "${LOG_DIR}/step2.HGDP_1KG.unrelated_reference_PCA.log" 2>&1

# --------------------------------------------------------------
# Step 3: Self-project all HGDP-1000G samples onto the PCA axes
#         defined by the unrelated reference subset.
#
# Use this .sscore file for plotting alongside CARTaGENE projection.
# --------------------------------------------------------------
# PC columns in .sscore start at column 6 (after IID, FID, and covariates)
PC_START_COL=6
PC_END_COL=$((5 + N_PCS))

echo "[$(date)] Step 3: Self-projecting all HGDP-1000G reference samples"

plink2 \
  --pfile "$REF_MERGED_PREFIX" \
  --read-freq "${REF_PCA_PREFIX}.acount" \
  --score "${REF_PCA_PREFIX}.eigenvec.allele" 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums "${PC_START_COL}-${PC_END_COL}" \
  --threads "$THREADS" \
  --out "$REF_SELF_PROJECT_PREFIX" \
  > "${LOG_DIR}/step3.HGDP_1KG.all_reference_self_projection.log" 2>&1

echo "[$(date)] HGDP/1000G reference PCA complete."
echo "KING retained sample list:"
echo "  ${KING_PREFIX}.king.cutoff.in.id"
echo "KING excluded sample list:"
echo "  ${KING_PREFIX}.king.cutoff.out.id"
echo "Unrelated-reference native PCA:"
echo "  ${REF_PCA_PREFIX}.eigenvec"
echo "Reference eigenvalues:"
echo "  ${REF_PCA_PREFIX}.eigenval"
echo "Reference allele weights:"
echo "  ${REF_PCA_PREFIX}.eigenvec.allele"
echo "Reference allele counts:"
echo "  ${REF_PCA_PREFIX}.acount"
echo "All-reference self-projected scores:"
echo "  ${REF_SELF_PROJECT_PREFIX}.sscore"