#!/bin/bash
# SLURM job array: genotype the ancestry panel's AIM sites in each tumor CRAM
# with scripts/genotype_aims.sh (one array task per CRAM).
#
# Usage (from anywhere; paths below are absolute):
#   sbatch --array=0 scripts/genotype_aims.slurm             # test: first CRAM only
#   N=$(ls $CRAM_GLOB | wc -l)                               # then: all CRAMs
#   sbatch --array=0-$((N - 1))%8 scripts/genotype_aims.slurm
#
# Override any path/module below with `sbatch --export=ALL,VAR=value ...`.
#
#SBATCH --job-name=aims_genotype
#SBATCH --partition=g64                  # small, I/O-bound job: g64 nodes (12 CPUs/64G) are plenty
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2                # mpileup/call are single-threaded; 1 extra for the pipe
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/data04/projects04/MarianaBoroni/lbbc_members/Michelle/variant_analysis/logs/%x_%A_%a.out
#SBATCH --error=/data04/projects04/MarianaBoroni/lbbc_members/Michelle/variant_analysis/logs/%x_%A_%a.err

set -euo pipefail

BASE=${BASE:-/data04/projects04/MarianaBoroni/lbbc_members/Michelle/variant_analysis}
REPO_DIR=${REPO_DIR:-$BASE/Variants_analysis}
REF=${REF:-/databases/sources/iGenomes/Homo_sapiens/GATK/GRCh38/Homo_sapiens_assembly38.fasta}
PANEL=${PANEL:-/data04/projects04/MarianaBoroni/lbbc_members/Michelle/AIMs/canonical_AIMs_panel_hg38.tsv}
CRAM_GLOB=${CRAM_GLOB:-/data02/workflows/USERS/mmparedese/BAMS/run*/preprocessing/recalibrated/*/*.recal.cram}
OUT=${OUT:-$BASE/aims_genotypes}
BCFTOOLS_MODULE=${BCFTOOLS_MODULE:-bcftools}
HTSLIB_MODULE=${HTSLIB_MODULE:-htslib}

module load "$BCFTOOLS_MODULE" "$HTSLIB_MODULE"
echo "[env] $(bcftools --version | head -1); $(bgzip --version | head -1)" >&2

mkdir -p "$OUT"
SITES=$OUT/aims_sites.tsv.gz

# Build the sites file once; flock keeps concurrent array tasks from racing.
(
  flock 9
  if [ ! -s "$SITES.tbi" ]; then
    "$REPO_DIR/scripts/genotype_aims.sh" sites "$PANEL" "$REF" "$SITES"
  fi
) 9> "$OUT/.sites.lock"

mapfile -t CRAMS < <(ls $CRAM_GLOB)
TASK=${SLURM_ARRAY_TASK_ID:-0}
[ "$TASK" -lt "${#CRAMS[@]}" ] || { echo "array index $TASK >= ${#CRAMS[@]} CRAMs" >&2; exit 1; }
CRAM=${CRAMS[$TASK]}
SAMPLE=$(basename "$CRAM" .recal.cram)
echo "[task $TASK/${#CRAMS[@]}] $SAMPLE <- $CRAM" >&2

"$REPO_DIR/scripts/genotype_aims.sh" call "$CRAM" "$REF" "$SITES" "$OUT/$SAMPLE.aims.vcf.gz"
