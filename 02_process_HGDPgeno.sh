#!/bin/bash

#SBATCH --account=rrg-vmooser
#SBATCH --job-name=02_process_HGDPgeno
#SBATCH --output=02_process_HGDPgeno.out
#SBATCH --error=02_process_HGDPgeno.err
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G

set -euo pipefail

# --------------------------------------------------------------
# Script efficiency (60912782)
# State: TIMEOUT (exit code 0)
# Nodes: 1
# Cores per node: 2
# CPU Utilized: 1-14:26:06
# CPU Efficiency: 48.96% of 3-06:30:02 core-walltime
# Job Wall-clock time: 1-15:15:01
# Memory Utilized: 1.18 GB
# Memory Efficiency: 9.83% of 12.00 GB

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
# Human Genome Diversity Project with 1000 Genomes supplement
IN_DIR="/lustre06/project/6061810/shared/HGDP_1KG/unfiltered_vcfs"

# One-column sample list (one HGDP-1000G sample ID per line)
# Includes only reference samples that passed QC (n = 4117)
REF_SAMPLES="/lustre06/project/6061810/chanalex/phd_thesis/Data/hgdp_1kg.final_keep_samples.txt"

# --------------------------------------------------------------
# Output directories
# --------------------------------------------------------------
BASE_OUT="/lustre07/scratch/chanalex/HGDP-1KG/PCA_QC"

REF_PASS_VCF_DIR="${BASE_OUT}/01_PASS_biallelic_sample_subset_vcf"
RAW_PGEN_DIR="${BASE_OUT}/02_raw_pgen"
SMISS_DIR="${BASE_OUT}/03_sample_missingness"
QC_PGEN_DIR="${BASE_OUT}/04_qc_pgen"
QC_VCF_DIR="${BASE_OUT}/05_qc_vcf_for_CARTaGENE_intersection"
LOG_DIR="${BASE_OUT}/logs"

mkdir -p \
  "$REF_PASS_VCF_DIR" \
  "$RAW_PGEN_DIR" \
  "$SMISS_DIR" \
  "$QC_PGEN_DIR" \
  "$QC_VCF_DIR" \
  "$LOG_DIR"

# Check that reference sample list exists and is not empty
if [[ ! -s "$REF_SAMPLES" ]]; then
  echo "ERROR: Reference sample list missing or empty: $REF_SAMPLES" >&2
  exit 1
fi

echo "Reference samples in keep list:"
wc -l "$REF_SAMPLES"

# --------------------------------------------------------------
# Parameters
# --------------------------------------------------------------
CHROMOSOMES=($(seq -f "chr%g" 1 22))

# Set thresholds
MAX_SAMPLE_MISSINGNESS=0.05   # 5% genotype missingness per sample (remove samples with call rate < 95%)
MAX_VARIANT_MISSINGNESS=0.05  # 5% genotype missingness per variant (remove variants missing in > 5% of samples)
MIN_MAF=0.01                  # Minor allele frequency > 1%
HWE_PVAL=1e-6                 # Variants which depart Hardy-Weinberg equilibrium (p-value < 1e-6)

# --------------------------------------------------------------
# Step 0: Confirm sample overlap with chr1 before running full pipeline
# --------------------------------------------------------------
CHR_TEST="chr1"
RAW_TEST_VCF="${IN_DIR}/gnomad.genomes.v3.1.2.hgdp_tgp.${CHR_TEST}.vcf.bgz"

echo "[$(date)] Step 0: Checking reference sample overlap on ${CHR_TEST}"

bcftools query -l "$RAW_TEST_VCF" \
  | sort \
  > "${BASE_OUT}/${CHR_TEST}.samples_in_vcf.txt"

sort "$REF_SAMPLES" \
  > "${BASE_OUT}/samples_requested.txt"

comm -12 \
  "${BASE_OUT}/samples_requested.txt" \
  "${BASE_OUT}/${CHR_TEST}.samples_in_vcf.txt" \
  > "${BASE_OUT}/samples_found_in_${CHR_TEST}.txt"

echo "Requested reference samples:"
wc -l "${BASE_OUT}/samples_requested.txt"

echo "Reference samples found in ${CHR_TEST} VCF:"
wc -l "${BASE_OUT}/samples_found_in_${CHR_TEST}.txt"

N_FOUND=$(wc -l < "${BASE_OUT}/samples_found_in_${CHR_TEST}.txt")

if [[ "$N_FOUND" -eq 0 ]]; then
  echo "ERROR: None of the requested reference samples were found in ${RAW_TEST_VCF}." >&2
  echo "Check whether REF_SAMPLES contains the correct sample IDs." >&2
  exit 1
fi

# Use the confirmed overlapping sample list for all downstream bcftools commands
REF_SAMPLES_CONFIRMED="${BASE_OUT}/samples_found_in_${CHR_TEST}.txt"

# --------------------------------------------------------------
# Step 1: Restrict raw HGDP-1000G VCFs to:
#         - those that passed hard QC filters in gnomAD
#         - FILTER=PASS variants (remove VQSR failure, AC0/low confidence, 
#           high inbreeding coefficient, etc.)
#         - biallelic SNPs only
#
#         - Keep genotype-only outputs for PCA
#
# Write BCF (instead of VCF) to save disk space and speed up processing
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do

  RAW_VCF="${IN_DIR}/gnomad.genomes.v3.1.2.hgdp_tgp.${CHR}.vcf.bgz"
  OUT_BCF="${REF_PASS_VCF_DIR}/${CHR}.HGDP_1KG.PASS.biallelic.samples.GTonly.bcf"
  LOG_FILE="${LOG_DIR}/${CHR}.step1.bcftools_PASS_biallelic_samples_GTonly.log"

  echo "[$(date)] Step 1: Processing $CHR with bcftools into GT-only BCF"

  # Check that input VCF exists and is indexed
  if [[ ! -f "$RAW_VCF" ]]; then
    echo "ERROR: Missing input VCF: $RAW_VCF" >&2
    exit 1
  fi

  # Check that VCF is indexed (either .tbi or .csi), if not, create an index
  if [[ ! -f "${RAW_VCF}.tbi" && ! -f "${RAW_VCF}.csi" ]]; then
    echo "Indexing raw VCF: $RAW_VCF"
    bcftools index -t "$RAW_VCF"
  fi

  # Check if output BCF already exists and is indexed; if so, skip processing
  if [[ -s "$OUT_BCF" && -s "${OUT_BCF}.csi" ]]; then
    echo "[$(date)] $CHR already completed; skipping."
    ### NOTE: for SLURM array jobs (running one chromosome per job), change "continue" with "exit 0"
    continue
  fi

  {
    bcftools view \
      --threads "$THREADS" \
      -S "$REF_SAMPLES_CONFIRMED" \
      -f PASS \
      -m2 -M2 \
      -v snps \
      -Ou \
      "$RAW_VCF" | \
    bcftools annotate \
      --threads "$THREADS" \
      -x INFO,FORMAT \
      -Ob \
      -o "$OUT_BCF" \
      -

    bcftools index "$OUT_BCF"
  } > "$LOG_FILE" 2>&1

done

# --------------------------------------------------------------
# Step 2: Convert filtered reference VCFs to PLINK2 PGEN
#         Use consistent variant IDs (i.e., CHR:POS:REF:ALT)
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do

  IN_BCF="${REF_PASS_VCF_DIR}/${CHR}.HGDP_1KG.PASS.biallelic.samples.GTonly.bcf"
  OUT_PREFIX="${RAW_PGEN_DIR}/${CHR}.HGDP_1KG.raw"
  LOG_FILE="${LOG_DIR}/${CHR}.step2.convert_to_pgen.log"

  echo "[$(date)] Step 2: Converting $CHR to PGEN"

  plink2 \
    --bcf "$IN_BCF" \
    --const-fid 0 \
    --snps-only just-acgt \
    --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' \
    --rm-dup exclude-all \
    --make-pgen \
    --threads "$THREADS" \
    --out "$OUT_PREFIX" \
    > "$LOG_FILE" 2>&1

done

# --------------------------------------------------------------
# Step 3: Compute per-sample missingness per chromosome
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do

  IN_PREFIX="${RAW_PGEN_DIR}/${CHR}.HGDP_1KG.raw"
  OUT_PREFIX="${SMISS_DIR}/${CHR}.HGDP_1KG"

  echo "[$(date)] Step 3: Computing reference sample missingness for $CHR"

  plink2 \
    --pfile "$IN_PREFIX" \
    --missing sample-only \
    --threads "$THREADS" \
    --out "$OUT_PREFIX" \
    > "${LOG_DIR}/${CHR}.step3.reference_smiss.log" 2>&1

done

# --------------------------------------------------------------
# Step 4: Aggregate reference sample missingness genome-wide
# --------------------------------------------------------------
GENOMEWIDE_REF_SMISS="${SMISS_DIR}/HGDP_1KG.genomewide.smiss"
KEEP_REF_SAMPLES="${SMISS_DIR}/HGDP_1KG.keep_samples.mind${MAX_SAMPLE_MISSINGNESS}.txt"
REMOVE_REF_SAMPLES="${SMISS_DIR}/HGDP_1KG.remove_samples.mind${MAX_SAMPLE_MISSINGNESS}.txt"

echo "[$(date)] Step 4: Aggregating reference genome-wide sample missingness"

awk '
BEGIN { OFS="\t" }

FNR == 1 {
  if ($1 == "#FID") {
    iid_col = 2
    miss_col = 3
    obs_col = 4
  } else if ($1 == "#IID") {
    iid_col = 1
    miss_col = 2
    obs_col = 3
  } else {
    print "ERROR: Unexpected .smiss header in " FILENAME ": " $0 > "/dev/stderr"
    exit 1
  }
  next
}

{
  iid = $iid_col
  miss[iid] += $miss_col
  obs[iid]  += $obs_col
}

END {
  print "#IID", "MISSING_CT", "OBS_CT", "F_MISS"
  for (iid in miss) {
    fmiss = (obs[iid] > 0 ? miss[iid] / obs[iid] : "NA")
    print iid, miss[iid], obs[iid], fmiss
  }
}
' "${SMISS_DIR}"/chr*.HGDP_1KG.smiss \
  | sort -k1,1 \
  > "$GENOMEWIDE_REF_SMISS"

awk -v max_miss="$MAX_SAMPLE_MISSINGNESS" '
NR > 1 && $4 <= max_miss { print $1 }
' "$GENOMEWIDE_REF_SMISS" > "$KEEP_REF_SAMPLES"

awk -v max_miss="$MAX_SAMPLE_MISSINGNESS" '
NR > 1 && $4 > max_miss { print $1 }
' "$GENOMEWIDE_REF_SMISS" > "$REMOVE_REF_SAMPLES"

echo "Reference samples kept after missingness:"
wc -l "$KEEP_REF_SAMPLES"

echo "Reference samples removed after missingness:"
wc -l "$REMOVE_REF_SAMPLES"

# --------------------------------------------------------------
# Step 5: Apply reference sample + variant QC
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do

  IN_PREFIX="${RAW_PGEN_DIR}/${CHR}.HGDP_1KG.raw"
  OUT_PREFIX="${QC_PGEN_DIR}/${CHR}.HGDP_1KG.QC"
  LOG_FILE="${LOG_DIR}/${CHR}.step5.reference_sample_variant_QC.log"

  echo "[$(date)] Step 5: Applying reference sample + variant QC for $CHR"

  plink2 \
    --pfile "$IN_PREFIX" \
    --keep "$KEEP_REF_SAMPLES" \
    --geno "$MAX_VARIANT_MISSINGNESS" \
    --maf "$MIN_MAF" \
    --hwe "$HWE_PVAL" 0.001 midp keep-fewhet \
    --make-pgen \
    --threads "$THREADS" \
    --out "$OUT_PREFIX" \
    > "$LOG_FILE" 2>&1

done

# --------------------------------------------------------------
# Step 6: Export QCed HGDP-1000G VCFs for intersection with
#         CARTaGENE QC VCFs
# --------------------------------------------------------------
for CHR in "${CHROMOSOMES[@]}"; do

  IN_PREFIX="${QC_PGEN_DIR}/${CHR}.HGDP_1KG.QC"
  OUT_PREFIX="${QC_VCF_DIR}/${CHR}.HGDP_1KG.QC"
  LOG_FILE="${LOG_DIR}/${CHR}.step6.export_reference_QC_vcf.log"

  echo "[$(date)] Step 6: Exporting QCed HGDP/1000G VCF for $CHR"

  plink2 \
    --pfile "$IN_PREFIX" \
    --export vcf bgz \
    --threads "$THREADS" \
    --out "$OUT_PREFIX" \
    > "$LOG_FILE" 2>&1

  bcftools index -t "${OUT_PREFIX}.vcf.gz"

done