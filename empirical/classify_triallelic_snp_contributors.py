#!/usr/bin/env python3

import sys
import subprocess
from collections import Counter

if len(sys.argv) != 2:
    sys.exit("Usage: python3 classify_triallelic_snp_contributors.py <vcf.gz>")

vcf = sys.argv[1]


def run_bcftools(vcf):
    view_cmd = [
        "bcftools", "view",
        "-m3", "-M3",
        "-i", "strlen(REF)=1 && strlen(ALT[0])=1 && strlen(ALT[1])=1",
        vcf,
    ]

    query_cmd = [
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n",
    ]

    p1 = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, text=True)
    p2 = subprocess.Popen(query_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
    p1.stdout.close()

    for line in p2.stdout:
        yield line.rstrip("\n")

    p2.stdout.close()


def gt_to_alleles(gt, alleles):
    gt = gt.replace("|", "/")

    if gt in {"./.", ".", ".|."}:
        return set()

    observed = set()

    for idx in gt.split("/"):
        if idx == ".":
            continue
        observed.add(alleles[int(idx)])

    return observed


def classify_two_sample_site(sample_alleles):
    """
    sample_alleles:
        dict: sample -> set(observed alleles)

    Categories are designed for two diploid samples, e.g.:

        REFIND = A/G
        S1     = A/T

    which gives an emergent triallelic SNP.
    """

    names = list(sample_alleles.keys())

    if len(names) != 2:
        return "not_two_samples"

    s1, s2 = names
    a1 = sample_alleles[s1]
    a2 = sample_alleles[s2]

    union = a1 | a2

    if len(union) < 3:
        return "triallelic_record_but_genotypes_show_less_than_3_alleles"

    shared = a1 & a2
    private1 = a1 - a2
    private2 = a2 - a1

    # Classic clean case:
    # REFIND = A/G
    # S1     = A/T
    # Both heterozygous, sharing one allele and each contributing one private allele.
    if len(a1) == 2 and len(a2) == 2 and len(shared) == 1 and len(private1) == 1 and len(private2) == 1:
        return "both_heterozygous_with_distinct_private_alleles"

    # One sample is homozygous for one allele;
    # the other is heterozygous for two different alleles.
    # Example:
    # REFIND = G/G
    # S1     = C/A
    if len(a1) == 1 and len(a2) == 2 and len(union) == 3:
        return f"{s1}_homozygous_one_allele_{s2}_heterozygous_two_other_alleles"

    if len(a2) == 1 and len(a1) == 2 and len(union) == 3:
        return f"{s2}_homozygous_one_allele_{s1}_heterozygous_two_other_alleles"

    # One sample has both alleles of the other plus an additional allele.
    # Rare for diploids, but useful if genotype parsing sees unusual calls.
    if len(private1) >= 1 and len(private2) == 0:
        return f"only_{s1}_contributes_private_allele"

    if len(private2) >= 1 and len(private1) == 0:
        return f"only_{s2}_contributes_private_allele"

    # General valid triallelic configuration.
    if len(union) == 3:
        return "other_valid_triallelic_configuration"

    return "other"


counts = Counter()
examples_by_category = {}

for line in run_bcftools(vcf):
    parts = line.split("\t")

    chrom, pos, ref, alts = parts[:4]
    sample_fields = parts[4:]

    alleles = [ref] + alts.split(",")

    sample_alleles = {}

    for field in sample_fields:
        sample, gt = field.split("=")
        sample_alleles[sample] = gt_to_alleles(gt, alleles)

    category = classify_two_sample_site(sample_alleles)
    counts[category] += 1

    if category not in examples_by_category:
        examples_by_category[category] = []

    if len(examples_by_category[category]) < 5:
        examples_by_category[category].append(
            (chrom, pos, ref, alts, sample_alleles)
        )


print("\nCategory counts")
print("---------------")
for category, count in counts.most_common():
    print(f"{category}\t{count}")


print("\nExamples by category")
print("--------------------")
for category, examples in examples_by_category.items():
    print(f"\n## {category}")
    for chrom, pos, ref, alts, sample_alleles in examples:
        print(f"{chrom}:{pos} {ref}>{alts} {sample_alleles}")


print("\nSummary")
print("-------")
total = sum(counts.values())
valid = total - counts.get("triallelic_record_but_genotypes_show_less_than_3_alleles", 0)

print(f"Total strict triallelic SNP records: {total}")
print(f"Valid genotype-supported triallelic sites: {valid}")

if total > 0:
    print(f"Fraction genotype-supported: {valid / total:.4f}")

core = counts.get("both_heterozygous_with_distinct_private_alleles", 0)
if total > 0:
    print(f"Fraction classic REFIND/focal distinct-private-allele pattern: {core / total:.4f}")
