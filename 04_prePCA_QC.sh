#!/bin/bash

#SBATCH --account=rrg-vmooser
#SBATCH --job-name=04_prePCA_QC
#SBATCH --output=04_prePCA_QC.out
#SBATCH --error=04_prePCA_QC.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

# --------------------------------------------------------------
# Script efficiency (61171399)
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
BASE_OUT="/lustre07/scratch/chanalex/CARTaGENE_HGDP-1KG/PCA_projection"

# Merged and QCed datasets from previous steps
REF_MERGED_PREFIX="${BASE_OUT}/09_merged_autosomes/HGDP_1KG.QC.shared_LDpruned.autosomes"
CAG_MERGED_PREFIX="${BASE_OUT}/09_merged_autosomes/CARTaGENE.QC.shared_LDpruned.autosomes"

# --------------------------------------------------------------
# Output directories
# --------------------------------------------------------------
QC_DIR="${BASE_OUT}/10_prePCA_sample_QC"
LOG_DIR="${BASE_OUT}/logs"

mkdir -p "$QC_DIR" "$LOG_DIR"

# --------------------------------------------------------------
# Step 1: Heterozygosity diagnostics
#
# Following: https://www.cog-genomics.org/plink/2.0/basic_stats#het
# PLINK2 --het computes observed/expected het counts and an F coefficient.
# PLINK notes this is best performed on a variant set in approximate
# HWE and linkage equilibrium.
# --------------------------------------------------------------
echo "[$(date)] Step 1: HGDP-1000G heterozygosity diagnostics"

plink2 \
  --pfile "$REF_MERGED_PREFIX" \
  --het \
  --threads "$THREADS" \
  --out "${QC_DIR}/HGDP_1KG.reference.het" \
  > "${LOG_DIR}/04.step1.HGDP_1KG.het.log" 2>&1

echo "[$(date)] Step 1: CARTaGENE heterozygosity diagnostics"

plink2 \
  --pfile "$CAG_MERGED_PREFIX" \
  --het \
  --threads "$THREADS" \
  --out "${QC_DIR}/CARTaGENE.projected_target.het" \
  > "${LOG_DIR}/04.step1.CARTaGENE.het.log" 2>&1

# --------------------------------------------------------------
# Step 2: Missingness diagnostics on the final SNP set
#
# This is only a check. Sample missingness filters already 
# applied upstream.
# --------------------------------------------------------------

echo "[$(date)] Step 2: HGDP-1000G missingness diagnostics"

plink2 \
  --pfile "$REF_MERGED_PREFIX" \
  --missing sample-only \
  --threads "$THREADS" \
  --out "${QC_DIR}/HGDP_1KG.reference.final_missingness" \
  > "${LOG_DIR}/04.step2.HGDP_1KG.final_missingness.log" 2>&1

echo "[$(date)] Step 2: CARTaGENE missingness diagnostics"

plink2 \
  --pfile "$CAG_MERGED_PREFIX" \
  --missing sample-only \
  --threads "$THREADS" \
  --out "${QC_DIR}/CARTaGENE.projected_target.final_missingness" \
  > "${LOG_DIR}/04.step2.CARTaGENE.final_missingness.log" 2>&1

# --------------------------------------------------------------
# Step 3: Optional global heterozygosity outlier lists
#
# Use for diagnostics, NOT automatic exclusions.
# For multi-ancestry data, global heterozygosity thresholds can flag
# ancestry differences rather than technical artefacts.
# --------------------------------------------------------------

for PREFIX in HGDP_1KG.reference CARTaGENE.projected_target; do

  HET_FILE="${QC_DIR}/${PREFIX}.het.het"
  OUTLIER_FILE="${QC_DIR}/${PREFIX}.het_outliers.global_3SD.txt"
  SUMMARY_FILE="${QC_DIR}/${PREFIX}.het_summary.txt"

  awk '
  BEGIN { OFS="\t" }
  NR == 1 {
    for (i=1; i<=NF; i++) {
      if ($i == "F") f_col=i
      if ($i == "#IID" || $i == "IID") iid_col=i
    }
    next
  }
  $f_col != "nan" && $f_col != "NA" {
    n++
    f[n]=$f_col
    iid[n]=$iid_col
    sum += $f_col
    sumsq += ($f_col * $f_col)
  }
  END {
    mean = sum / n
    sd = sqrt((sumsq - (sum * sum / n)) / (n - 1))
    lower = mean - 3 * sd
    upper = mean + 3 * sd

    print "N", n
    print "Mean_F", mean
    print "SD_F", sd
    print "Lower_3SD", lower
    print "Upper_3SD", upper

    for (i=1; i<=n; i++) {
      if (f[i] < lower || f[i] > upper) {
        print iid[i], f[i] > "'$OUTLIER_FILE'"
      }
    }
  }
  ' "$HET_FILE" > "$SUMMARY_FILE"

done

echo "[$(date)] Pre-PCA sample QC diagnostics complete."
echo "Outputs:"
echo "  $QC_DIR"