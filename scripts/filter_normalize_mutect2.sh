#!/usr/bin/env bash
# Filter a Mutect2+VEP VCF to PASS/germline/panel_of_normals/combinations and
# normalize it (split multiallelics, left-align against the reference) before
# it is fed to `tumoronly run`. See docs/GETTING_STARTED.md and
# docs/FILTERING_STRATEGY.md ("Caller FILTER as graded evidence") for why
# these specific FILTER categories are kept, and why normalization happens
# here rather than inside tumoronly.
#
# Usage:
#   scripts/filter_normalize_mutect2.sh <input.vcf.gz> <reference.fasta> <output_prefix>
#
# Produces:
#   <output_prefix>.PASS.vcf.gz          retained (PASS/germline/panel_of_normals/combinations), normalized
#   <output_prefix>.PASS.vcf.gz.tbi
#   <output_prefix>.nonPASS.vcf.gz       everything else, kept (not deleted), NOT normalized
#
# Requires bcftools (and bgzip/tabix, bundled with it) on PATH.

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <input.vcf.gz> <reference.fasta> <output_prefix>" >&2
  exit 1
fi

input_vcf="$1"
fasta="$2"
out_prefix="$3"

command -v bcftools >/dev/null 2>&1 || {
  echo "bcftools not found on PATH. Activate an environment that provides it, e.g.:" >&2
  echo "  export PATH=/home/mmparedese/miniconda3/envs/vcfman/bin:\$PATH" >&2
  exit 1
}
[ -f "$input_vcf" ] || { echo "Input VCF not found: $input_vcf" >&2; exit 1; }
[ -f "$fasta" ] || { echo "Reference FASTA not found: $fasta" >&2; exit 1; }

accepted_filter='FILTER="PASS" || FILTER="germline" || FILTER="panel_of_normals" || FILTER="germline;panel_of_normals"'

echo "[1/3] Splitting by FILTER (PASS/germline/panel_of_normals/combinations vs. everything else)..." >&2
bcftools view -i "$accepted_filter" "$input_vcf" -Oz -o "${out_prefix}.PASS.unnormalized.vcf.gz"
bcftools view -e "$accepted_filter" "$input_vcf" -Oz -o "${out_prefix}.nonPASS.vcf.gz"
bcftools index -t "${out_prefix}.nonPASS.vcf.gz"

echo "[2/3] Normalizing retained variants (split multiallelics, left-align vs. reference)..." >&2
bcftools norm -f "$fasta" -m -any "${out_prefix}.PASS.unnormalized.vcf.gz" -Oz -o "${out_prefix}.PASS.vcf.gz"
bcftools index -t "${out_prefix}.PASS.vcf.gz"
rm -f "${out_prefix}.PASS.unnormalized.vcf.gz"

echo "[3/3] Done." >&2
echo "  retained (normalized): ${out_prefix}.PASS.vcf.gz    ($(bcftools view -H "${out_prefix}.PASS.vcf.gz" | wc -l) variants)" >&2
echo "  excluded (kept, raw):  ${out_prefix}.nonPASS.vcf.gz ($(bcftools view -H "${out_prefix}.nonPASS.vcf.gz" | wc -l) variants)" >&2
