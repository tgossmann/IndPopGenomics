#!/usr/bin/env bash
set -euo pipefail

# Paths to the reference assembly and the mapped reads are set in config.sh
source "$(dirname "$0")/config.sh"

OUT="realdata_freebayes_highcov_${CHR}"
mkdir -p "$OUT"/{vcf,stats,logs,triallelic}

if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

call_freebayes () {
    NAME="$1"
    shift
    BAMS=("$@")

    echo "Calling $NAME"

    freebayes \
        -f "$REF" \
        -r "$CHR" \
        --genotype-qualities \
        --min-base-quality 20 \
        --min-mapping-quality 30 \
        --min-coverage 8 \
        --use-best-n-alleles 6 \
        "${BAMS[@]}" \
        > "$OUT/vcf/${NAME}.raw.vcf" \
        2> "$OUT/logs/${NAME}.freebayes.log"

    bgzip -f "$OUT/vcf/${NAME}.raw.vcf"
    tabix -f -p vcf "$OUT/vcf/${NAME}.raw.vcf.gz"

    # Keep SNP records, including multiallelic SNP records.
    bcftools view \
        -v snps \
        "$OUT/vcf/${NAME}.raw.vcf.gz" \
        -Oz -o "$OUT/vcf/${NAME}.snps.multiallelic.vcf.gz"

    tabix -f -p vcf "$OUT/vcf/${NAME}.snps.multiallelic.vcf.gz"

    # Hard filter.
    bcftools filter \
        -e 'QUAL<30 || INFO/DP<8' \
        "$OUT/vcf/${NAME}.snps.multiallelic.vcf.gz" \
        -Oz -o "$OUT/vcf/${NAME}.snps.multiallelic.filtered.vcf.gz"

    tabix -f -p vcf "$OUT/vcf/${NAME}.snps.multiallelic.filtered.vcf.gz"

    # Strict single-base triallelic SNPs only:
    # REF length 1, exactly two ALT alleles, both ALT length 1.
    bcftools view \
        -m3 -M3 \
        -i 'strlen(REF)=1 && strlen(ALT[0])=1 && strlen(ALT[1])=1' \
        "$OUT/vcf/${NAME}.snps.multiallelic.filtered.vcf.gz" \
        -Oz -o "$OUT/triallelic/${NAME}.strict_triallelic_snps.vcf.gz"

    tabix -f -p vcf "$OUT/triallelic/${NAME}.strict_triallelic_snps.vcf.gz"

    # Position list
    bcftools query \
        -f '%CHROM\t%POS\n' \
        "$OUT/triallelic/${NAME}.strict_triallelic_snps.vcf.gz" \
        | sort -k1,1 -k2,2n \
        > "$OUT/triallelic/${NAME}.strict_triallelic_snps.pos"

    # Count
    N=$(bcftools view -H "$OUT/triallelic/${NAME}.strict_triallelic_snps.vcf.gz" | wc -l | tr -d ' ')
    echo -e "${NAME}\t${N}" >> "$OUT/triallelic/strict_triallelic_counts.tsv"
}

rm -f "$OUT/triallelic/strict_triallelic_counts.tsv"
echo -e "analysis\tstrict_triallelic_snps" > "$OUT/triallelic/strict_triallelic_counts.tsv"

# -----------------------------
# A. High-coverage focal individuals only
# -----------------------------
call_freebayes "H1_only" "$H1"
call_freebayes "H2_only" "$H2"
call_freebayes "H3_only" "$H3"
call_freebayes "H4_only" "$H4"

# -----------------------------
# B. Reference individual only
# -----------------------------
call_freebayes "REFIND_only" "$REFIND"

# -----------------------------
# C. Reference individual + one high-coverage focal
# -----------------------------
call_freebayes "REFIND_plus_H1" "$REFIND" "$H1"
call_freebayes "REFIND_plus_H2" "$REFIND" "$H2"
call_freebayes "REFIND_plus_H3" "$REFIND" "$H3"
call_freebayes "REFIND_plus_H4" "$REFIND" "$H4"

# -----------------------------
# D. All four high-coverage focals without reference donor
# -----------------------------
call_freebayes "ALL4_highcov_no_REFIND" "$H1" "$H2" "$H3" "$H4"

# -----------------------------
# E. Reference donor + all four high-coverage focals
# -----------------------------
call_freebayes "REFIND_plus_ALL4_highcov" "$REFIND" "$H1" "$H2" "$H3" "$H4"

echo "Done."
echo "Output folder: $OUT"
echo "Triallelic counts:"
cat "$OUT/triallelic/strict_triallelic_counts.tsv"
