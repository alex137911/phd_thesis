#!/bin/bash

#SBATCH --account=rrg-vmooser
#SBATCH --job-name=01_process_CAGgeno
#SBATCH --output=01_process_CAGgeno.out
#SBATCH --error=01_process_CAGgeno.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

# Exit on error, treat unset variables as errors, and prevent errors in 
# pipelines from being masked
set -euo pipefail

# --------------------------------------------------------------
# Script efficiency (60671816)
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 2
# CPU Utilized: 04:20:42
# CPU Efficiency: 52.74% of 08:14:18 core-walltime
# Job Wall-clock time: 04:07:09
# Memory Utilized: 28.44 MB
# Memory Efficiency: 0.35% of 8.00 GB

# --------------------------------------------------------------
# Load necessary modules
# --------------------------------------------------------------
module purge
module load StdEnv/2023
module load gcc/12.3
module load bcftools/1.19
module load plink/2.00a5.8

# Set threads for parallel processing
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# --------------------------------------------------------------
# Directories
# --------------------------------------------------------------
# Input directory containing raw VCF files
IN_DIR="/lustre06/project/6061810/CERC_Private/Geno/CARTaGENE/Array/Pre-imputation"

# Output directories for each step of the pipeline
BASE_OUT="/lustre07/scratch/chanalex/CARTaGENE/PCA_QC"

RAW_PGEN_DIR="${BASE_OUT}/01_biallelic_raw_pgen"
SMISS_DIR="${BASE_OUT}/02_sample_missingness"
QC_PGEN_DIR="${BASE_OUT}/03_qc_pgen"
QC_VCF_DIR="${BASE_OUT}/04_qc_vcf_for_reference_intersection"
LOG_DIR="${BASE_OUT}/logs"

mkdir -p \
  "$RAW_PGEN_DIR" \
  "$SMISS_DIR" \
  "$QC_PGEN_DIR" \
  "$QC_VCF_DIR" \
  "$LOG_DIR"

# --------------------------------------------------------------
# Parameters & Thresholds
# --------------------------------------------------------------
# Chromosomes to process (autosomes only for PCA)
CHROMOSOMES=($(seq -f "chr%g" 1 22))

# Set thresholds
MAX_SAMPLE_MISSINGNESS=0.05   # 5% genotype missingness per sample (remove samples with call rate < 95%)
MAX_VARIANT_MISSINGNESS=0.05  # 5% genotype missingness per variant (remove variants missing in > 5% of samples)
MIN_MAF=0.01                  # Minor allele frequency > 1%
# HWE_PVAL=1e-6               # Variants which depart Hardy-Weinberg equilibrium (p-value < 1e-6)

# --------------------------------------------------------------
# Step 1: Convert each CARTaGENE chromosome .VCF to PLINK2 PGEN
#         Keep autosomal biallelic A/C/G/T SNPs only.
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do
  INPUT_VCF="${IN_DIR}/${CHR}.CARTaGENEv1.1.array.vcf.gz"
  OUT_PREFIX="${RAW_PGEN_DIR}/${CHR}.CARTaGENE.biallelic_raw"

  if [[ ! -f "$INPUT_VCF" ]]; then
    echo "ERROR: Missing input VCF: $INPUT_VCF" >&2
    exit 1
  fi

  echo "[$(date)] Step 1: Converting $CHR to biallelic SNP PGEN"

# PLINK2 command to convert VCF to PGEN, keeping only biallelic A/C/G/T SNPs
  plink2 \
    --vcf "$INPUT_VCF" \
    --const-fid 0 \
    --snps-only just-acgt \
    --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' \
    --rm-dup exclude-all \
    --make-pgen \
    --threads "$THREADS" \
    --out "$OUT_PREFIX" \
    > "${LOG_DIR}/${CHR}.step1.convert.log" 2>&1
done

# --------------------------------------------------------------
# Step 2: Compute per-sample missingness per chromosome
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do
  IN_PREFIX="${RAW_PGEN_DIR}/${CHR}.CARTaGENE.biallelic_raw"
  OUT_PREFIX="${SMISS_DIR}/${CHR}.CARTaGENE"

  echo "[$(date)] Step 2: Computing sample missingness for $CHR"

  plink2 \
    --pfile "$IN_PREFIX" \
    --missing sample-only \
    --threads "$THREADS" \
    --out "$OUT_PREFIX" \
    > "${LOG_DIR}/${CHR}.step2.smiss.log" 2>&1
done

# --------------------------------------------------------------
# Step 3: Aggregate sample missingness genome-wide
#         - creates one sample list used for all chromosomes
# --------------------------------------------------------------
GENOMEWIDE_SMISS="${SMISS_DIR}/CARTaGENE.genomewide.smiss"
KEEP_SAMPLES="${SMISS_DIR}/CARTaGENE.keep_samples.mind${MAX_SAMPLE_MISSINGNESS}.txt"
REMOVE_SAMPLES="${SMISS_DIR}/CARTaGENE.remove_samples.mind${MAX_SAMPLE_MISSINGNESS}.txt"

echo "[$(date)] Step 3: Aggregating genome-wide sample missingness"

awk '
BEGIN { OFS="\t" }
FNR == 1 { next }
{
  key = $1 "\t" $2
  miss[key] += $3
  obs[key]  += $4
}
END {
  print "#FID", "IID", "MISSING_CT", "OBS_CT", "F_MISS"
  for (key in miss) {
    split(key, a, "\t")
    fmiss = (obs[key] > 0 ? miss[key] / obs[key] : "NA")
    print a[1], a[2], miss[key], obs[key], fmiss
  }
}
' "${SMISS_DIR}"/chr*.CARTaGENE.smiss \
  | sort -k2,2 \
  > "$GENOMEWIDE_SMISS"

awk -v max_miss="$MAX_SAMPLE_MISSINGNESS" '
BEGIN { OFS="\t" }
NR > 1 && $5 <= max_miss { print $1, $2 }
' "$GENOMEWIDE_SMISS" > "$KEEP_SAMPLES"

awk -v max_miss="$MAX_SAMPLE_MISSINGNESS" '
BEGIN { OFS="\t" }
NR > 1 && $5 > max_miss { print $1, $2 }
' "$GENOMEWIDE_SMISS" > "$REMOVE_SAMPLES"

echo "Samples kept:"
wc -l "$KEEP_SAMPLES"

echo "Samples removed:"
wc -l "$REMOVE_SAMPLES"

# --------------------------------------------------------------
# Step 4: Apply one genome-wide sample keep list to all chromosomes
#         and filter variants by missingness and MAF
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do
  IN_PREFIX="${RAW_PGEN_DIR}/${CHR}.CARTaGENE.biallelic_raw"
  OUT_PREFIX="${QC_PGEN_DIR}/${CHR}.CARTaGENE.QC"

  echo "[$(date)] Step 4: Applying sample + variant QC for $CHR"

  plink2 \
    --pfile "$IN_PREFIX" \
    --keep "$KEEP_SAMPLES" \
    --geno "$MAX_VARIANT_MISSINGNESS" \
    --maf "$MIN_CAG_MAF" \
    --make-pgen \
    --threads "$THREADS" \
    --out "$OUT_PREFIX" \
    > "${LOG_DIR}/${CHR}.step4.CARTaGENE_QC.log" 2>&1
done

# --------------------------------------------------------------
# Step 5: Export QCed CARTaGENE VCFs for HGDP/1000G intersection.
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do
  IN_PREFIX="${QC_PGEN_DIR}/${CHR}.CARTaGENE.QC"
  OUT_PREFIX="${QC_VCF_DIR}/${CHR}.CARTaGENE.QC"

  echo "[$(date)] Step 5: Exporting QCed VCF for $CHR"

  plink2 \
    --pfile "$IN_PREFIX" \
    --export vcf bgz \
    --threads "$THREADS" \
    --out "$OUT_PREFIX" \
    > "${LOG_DIR}/${CHR}.step5.export_QC_vcf.log" 2>&1

  bcftools index -t "${OUT_PREFIX}.vcf.gz"
done

echo "[$(date)] CARTaGENE QC complete."
echo "QC PGEN files:"
echo "  $QC_PGEN_DIR"
echo "QC VCFs for HGDP/1000G intersection:"
echo "  $QC_VCF_DIR"