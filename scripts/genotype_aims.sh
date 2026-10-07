#!/usr/bin/env bash
# Genotype the ancestry-informative markers (AIMs) of the ancestry panel
# directly from an aligned CRAM/BAM, emitting EVERY covered panel site --
# including homozygous-reference (0/0) ones -- with per-allele read counts.
#
# Why: caller VCFs (Mutect2, sarek's `bcftools call -mv`) only contain sites
# with an ALT allele, so 0/0 AIMs are silently absent and the NNLS ancestry
# estimate is biased toward populations where the ALT allele is common. See
# docs/ANCESTRY.md. These per-sample VCFs are a separate ancestry track: they
# are NOT annotated (no VEP needed) and NOT passed through somatic filters.
#
# Usage:
#   scripts/genotype_aims.sh sites <aims_panel.tsv> <reference.fasta> <out_sites.tsv.gz>
#       Build the sites file once: biallelic SNVs from the panel as
#       CHROM<TAB>POS<TAB>REF,ALT, contig names matched to the FASTA
#       ("chr1" vs "1"), sorted in FASTA order, bgzipped + tabix-indexed.
#       The panel needs CHROM/POS/REF/ALT columns (case-insensitive; also
#       chromosome/position) or a canonical_key column (BUILD|CHROM|POS|REF|ALT).
#
#   scripts/genotype_aims.sh call <sample.cram|bam> <reference.fasta> <sites.tsv.gz> <out.vcf.gz>
#       Genotype one sample at the panel sites (bcftools mpileup | call -m
#       -C alleles, no -v), normalize against the reference, and index.
#       Sites with no read coverage are absent (= missing, not 0/0).
#
# Requires bcftools (with bgzip/tabix) on PATH.

set -euo pipefail

die() { echo "$*" >&2; exit 1; }

command -v bcftools >/dev/null 2>&1 || {
  echo "bcftools not found on PATH. Activate an environment that provides it, e.g.:" >&2
  echo "  export PATH=/home/mmparedese/miniconda3/envs/vcfman/bin:\$PATH" >&2
  exit 1
}

mode="${1:-}"
case "$mode" in
  sites)
    [ "$#" -eq 4 ] || die "Usage: $0 sites <aims_panel.tsv> <reference.fasta> <out_sites.tsv.gz>"
    panel="$2"; fasta="$3"; out="$4"
    [ -f "$panel" ] || die "Panel not found: $panel"
    [ -f "$fasta.fai" ] || die "FASTA index not found: $fasta.fai (run: samtools faidx $fasta)"
    case "$out" in *.tsv.gz) ;; *) die "Output must end in .tsv.gz: $out" ;; esac

    # contigs in FASTA order; whether the FASTA uses the "chr" prefix
    fai_has_chr=$(awk 'NR == 1 { print ($1 ~ /^chr/) ? 1 : 0 }' "$fasta.fai")

    awk -F'\t' -v OFS='\t' -v has_chr="$fai_has_chr" '
      FNR == NR { order[$1] = FNR; next }          # first file: the .fai
      FNR == 1 {
        for (i = 1; i <= NF; i++) {
          h = tolower($i); gsub(/^#/, "", h); gsub(/\r/, "", h)
          if (h == "chrom" || h == "chromosome" || h == "chr") c = i
          else if (h == "pos" || h == "position") p = i
          else if (h == "ref") r = i
          else if (h == "alt") a = i
          else if (h == "canonical_key") k = i
        }
        if (!(c && p && r && a) && !k) {
          print "panel header needs CHROM/POS/REF/ALT or canonical_key columns" > "/dev/stderr"; exit 2
        }
        next
      }
      {
        gsub(/\r/, "")
        if (c && p && r && a) { chrom = $c; pos = $p; ref = toupper($r); alt = toupper($a) }
        else { split($k, f, "|"); chrom = f[2]; pos = f[3]; ref = toupper(f[4]); alt = toupper(f[5]) }
        if (ref !~ /^[ACGT]$/ || alt !~ /^[ACGT]$/ || pos !~ /^[0-9]+$/) { skipped++; next }
        sub(/^chr/, "", chrom); if (chrom == "M") chrom = "MT"
        if (has_chr) chrom = (chrom == "MT") ? "chrM" : "chr" chrom
        if (!(chrom in order)) { missing++; next }
        key = chrom SUBSEP pos SUBSEP ref SUBSEP alt
        if (key in seen) next
        seen[key] = 1
        print order[chrom], chrom, pos, ref "," alt
        kept++
      }
      END {
        printf "[sites] kept=%d skipped_non_snv=%d contig_not_in_fasta=%d\n", kept, skipped, missing > "/dev/stderr"
      }' "$fasta.fai" "$panel" |
      sort -t$'\t' -k1,1n -k3,3n | cut -f2- | bgzip -c > "$out"
    tabix -f -s1 -b2 -e2 "$out"
    echo "[sites] written: $out ($(zcat "$out" | wc -l) sites)" >&2
    ;;

  call)
    [ "$#" -eq 5 ] || die "Usage: $0 call <sample.cram|bam> <reference.fasta> <sites.tsv.gz> <out.vcf.gz>"
    aln="$2"; fasta="$3"; sites="$4"; out="$5"
    [ -f "$aln" ] || die "Alignment not found: $aln"
    [ -f "$fasta" ] || die "Reference FASTA not found: $fasta"
    [ -f "$sites" ] && [ -f "$sites.tbi" ] || die "Sites file or its .tbi not found: $sites (run the 'sites' mode first)"
    case "$out" in *.vcf.gz) ;; *) die "Output must end in .vcf.gz: $out" ;; esac

    # -R: jump to panel sites via the index; -C alleles -T: force the panel's
    # REF/ALT; no -v, so 0/0 sites are kept; -A keeps the panel ALT on 0/0
    # sites (otherwise ALT="." and the site can't be matched to the panel by
    # canonical key). MAPQ/BASEQ >= 20 as in GATK-style QC.
    bcftools mpileup -f "$fasta" -R "$sites" -a FORMAT/AD,FORMAT/DP \
        -q 20 -Q 20 --max-depth 1000 -Ou "$aln" |
      bcftools call -m -A -C alleles -T "$sites" -Ou |
      bcftools norm -f "$fasta" -c ws -Oz -o "$out"
    bcftools index -t -f "$out"

    n_total=$(bcftools view -H "$out" | wc -l)
    n_homref=$(bcftools view -H -i 'GT="0/0"' "$out" | wc -l)
    echo "[call] written: $out (sites=$n_total, 0/0=$n_homref)" >&2
    ;;

  *)
    die "Usage: $0 sites <aims_panel.tsv> <reference.fasta> <out_sites.tsv.gz>
       $0 call  <sample.cram|bam> <reference.fasta> <sites.tsv.gz> <out.vcf.gz>"
    ;;
esac
