#!/usr/bin/env bash
set -euo pipefail

CHR="NC_031792.1"

OUT="realdata_freebayes_highcov_${CHR}"
STATDIR="$OUT/stats"
mkdir -p "$STATDIR"

# These should already exist from the high-coverage FreeBayes run
VCFDIR="$OUT/vcf"

for NAME in \
H1_only \
H2_only \
H3_only \
H4_only \
REFIND_only \
REFIND_plus_H1 \
REFIND_plus_H2 \
REFIND_plus_H3 \
REFIND_plus_H4 \
ALL4_highcov_no_REFIND \
REFIND_plus_ALL4_highcov
do
    echo "Processing $NAME"

    VCF="$VCFDIR/${NAME}.snps.multiallelic.filtered.vcf.gz"

    # Make biallelic branch for standard pi/TajimaD
    bcftools view \
        -m2 -M2 -v snps \
        "$VCF" \
        -Oz -o "$VCFDIR/${NAME}.snps.biallelic.filtered.vcf.gz"

    tabix -f -p vcf "$VCFDIR/${NAME}.snps.biallelic.filtered.vcf.gz"

    vcftools \
        --gzvcf "$VCFDIR/${NAME}.snps.biallelic.filtered.vcf.gz" \
        --window-pi 50000 \
        --window-pi-step 50000 \
        --out "$STATDIR/${NAME}"

    vcftools \
        --gzvcf "$VCFDIR/${NAME}.snps.biallelic.filtered.vcf.gz" \
        --TajimaD 50000 \
        --out "$STATDIR/${NAME}"
done

echo "Done. Stats written to $STATDIR"
