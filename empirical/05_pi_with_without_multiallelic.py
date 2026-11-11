#!/usr/bin/env python3
"""Revision 1 (R2.3): effect of including or excluding multiallelic SNP records on
nucleotide diversity, genome-wide (chromosome-wide) and in 50-kb windows.

Treatments, applied to <NAME>.snps.multiallelic.filtered.vcf.gz:
  biallelic  : records with exactly two alleles (REF + 1 ALT)  [as in the manuscript]
  + ref-induced : additionally records with more than two alleles at which the
                  sampled individuals carry only two alleles (e.g. reference A,
                  focal C/T); within the sample these are ordinary biallelic sites
  all        : additionally records at which the sampled individuals carry
               three or more alleles (genuinely multiallelic)
Per-site pi = n/(n-1) * (1 - sum p_i^2) over called haplotypes; window pi = sum / 50,000.
Only single-base SNP alleles are used.

Usage: python3 triallelic_pi_rev1.py <OUTDIR>
"""
import sys, gzip
from pathlib import Path
import numpy as np, pandas as pd
OUTDIR = Path(sys.argv[1]); W = 50_000
NAMES = ["H1_only","H2_only","H3_only","H4_only","REFIND_only",
         "REFIND_plus_H1","REFIND_plus_H2","REFIND_plus_H3","REFIND_plus_H4",
         "ALL4_highcov_no_REFIND","REFIND_plus_ALL4_highcov"]

def site_pi(alleles):
    n = len(alleles)
    if n < 2: return 0.0, 0
    _, c = np.unique(alleles, return_counts=True)
    p = c / n
    return n / (n - 1) * (1 - np.sum(p * p)), len(c)

rows, wins = [], []
for name in NAMES:
    acc = {}  # window -> [bi, refind, multi]
    counts = dict(bi=0, refinduced=0, multi=0)
    with gzip.open(OUTDIR / "vcf" / f"{name}.snps.multiallelic.filtered.vcf.gz", "rt") as fh:
        for line in fh:
            if line[0] == "#": continue
            f = line.rstrip("\n").split("\t")
            ref, alts = f[3], f[4].split(",")
            if len(ref) != 1 or any(len(a) != 1 for a in alts): continue
            al = []
            for g in f[9:]:
                gt = g.split(":")[0].replace("|", "/")
                al += [int(x) for x in gt.split("/") if x != "."]
            pi, k = site_pi(np.array(al))
            if k < 2: continue                      # not segregating among the samples
            b = (int(f[1]) - 1) // W * W + 1
            a = acc.setdefault(b, [0.0, 0.0, 0.0])
            if len(alts) == 1:
                a[0] += pi; counts["bi"] += 1
            elif k == 2:
                a[1] += pi; counts["refinduced"] += 1
            else:
                a[2] += pi; counts["multi"] += 1
    for b, (x, y, z) in acc.items():
        wins.append(dict(analysis=name, BIN_START=b, pi_biallelic=x / W,
                         pi_plus_refinduced=(x + y) / W, pi_all=(x + y + z) / W))
    rows.append(dict(analysis=name, sites_biallelic=counts["bi"],
                     sites_refinduced=counts["refinduced"], sites_multiallelic=counts["multi"]))

wdf = pd.DataFrame(wins)
cnt = pd.DataFrame(rows)
# restrict to the 137 windows used in Figure 4 (present in every vcftools output)
vt = pd.read_csv(OUTDIR / "stats" / "REFIND_plus_ALL4_highcov.windowed.pi", sep="\t")
common = set(vt["BIN_START"])
for n in NAMES:
    common &= set(pd.read_csv(OUTDIR / "stats" / f"{n}.windowed.pi", sep="\t")["BIN_START"])
wdf = wdf[wdf["BIN_START"].isin(common)]
summ = wdf.groupby("analysis")[["pi_biallelic","pi_plus_refinduced","pi_all"]].mean().reset_index()
summ = cnt.merge(summ, on="analysis")
summ["pct_change_refinduced"] = 100 * (summ.pi_plus_refinduced / summ.pi_biallelic - 1)
summ["pct_change_all"] = 100 * (summ.pi_all / summ.pi_biallelic - 1)
# check against vcftools biallelic pi
summ["pi_vcftools"] = [pd.read_csv(OUTDIR / "stats" / f"{n}.windowed.pi", sep="\t")
                       .query("BIN_START in @common")["PI"].mean() for n in summ.analysis]
summ = summ.set_index("analysis").loc[NAMES].reset_index()
outd = OUTDIR / "triallelic"; outd.mkdir(exist_ok=True)
wdf.to_csv(outd / "pi_with_without_multiallelic_windows_rev1.csv", index=False)
summ.to_csv(outd / "pi_with_without_multiallelic_summary_rev1.csv", index=False)
pd.set_option("display.width", 250)
print(summ.round(6).to_string(index=False))
