#!/bin/bash
# SLURM wrapper for scripts/run_ancestry_only.R: ancestry module only, on all
# VCFs matching INPUT_GLOB, written to OUTPUT. One job processes the VCFs one
# at a time (memory is bounded by the largest single VCF).
#
# Usage (set INPUT_GLOB/OUTPUT per test):
#   sbatch --export=ALL,INPUT_GLOB='/path/*/*.vcf.gz',OUTPUT=/path/out scripts/run_ancestry_only.slurm
#
#SBATCH --job-name=ancestry_only
#SBATCH --partition=g128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8                # data.table threads for VCF parsing
#SBATCH --mem=64G                        # one caller VCF (~5M records) at a time; AIM-only VCFs need <2G
#SBATCH --time=08:00:00
#SBATCH --output=/data04/projects04/MarianaBoroni/lbbc_members/Michelle/variant_analysis/logs/%x_%j.out
#SBATCH --error=/data04/projects04/MarianaBoroni/lbbc_members/Michelle/variant_analysis/logs/%x_%j.err

set -euo pipefail

BASE=${BASE:-/data04/projects04/MarianaBoroni/lbbc_members/Michelle/variant_analysis}
REPO_DIR=${REPO_DIR:-$BASE/Variants_analysis}
CONFIG=${CONFIG:-$BASE/example_ov_wgs.yml}
: "${INPUT_GLOB:?set INPUT_GLOB, e.g. --export=ALL,INPUT_GLOB='/path/*/*.vcf.gz'}"
: "${OUTPUT:?set OUTPUT, e.g. --export=ALL,OUTPUT=$BASE/results/ancestry_test}"

module load R/4.3.0
export R_LIBS_USER="/data04/tools/R/Rpackages/4.4.1"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_DATATABLE_NUM_THREADS=$SLURM_CPUS_PER_TASK

mapfile -t VCFS < <(ls $INPUT_GLOB)
[ "${#VCFS[@]}" -gt 0 ] || { echo "no VCFs match: $INPUT_GLOB" >&2; exit 1; }
echo "[ancestry_only] ${#VCFS[@]} VCF(s) -> $OUTPUT" >&2
printf '  %s\n' "${VCFS[@]}" >&2

# config paths (db/, config/) are relative to the repo, as in run_ov_wgs.slurm
cd "$REPO_DIR"
stdbuf -oL -eL Rscript scripts/run_ancestry_only.R --config "$CONFIG" --output "$OUTPUT" "${VCFS[@]}"
