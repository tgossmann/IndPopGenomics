#!/usr/bin/env python3
"""Revision 1, Figure 4: windowed pi, Watterson's theta and Tajima's D for one
focal individual analysed alone (2 haplotypes), against the pseudohaploid
reference (3 sequences), with the genotyped reference donor (4 haplotypes),
and for four focal individuals without (8) and with (10) the donor.

Usage: python3 plot_publication_popgen_highcov_rev1.py <OUTDIR> <H1|H2|H3|H4> <SRR_ID>

The pseudohaploid case is computed from <SAMPLE>_only.snps.biallelic.filtered.vcf.gz:
the reference base is treated as a third sequence, so every site at which the
focal individual carries a non-reference allele (0/1 or 1/1) is segregating among
the three sequences. For three sequences pi = theta_W = (2/3) * S / L at every
site, and Tajima's D is undefined.
"""
import sys, gzip
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

if len(sys.argv) != 4:
    sys.exit(__doc__)
OUTDIR, SAMPLE, SRR_ID = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
STATDIR, VCFDIR = OUTDIR / "stats", OUTDIR / "vcf"
PLOTDIR = OUTDIR / "plots_publication"; PLOTDIR.mkdir(exist_ok=True)
W = 50_000

CONFIGS = [  # key, analysis name, haplotypes, label
    ("focal",     f"{SAMPLE}_only",             2,  "Focal\nalone (2)"),
    ("pseudo",    None,                         3,  "Pseudo-\nhaploid (3)"),
    ("focal_ref", f"REFIND_plus_{SAMPLE}",      4,  "Focal +\ndonor (4)"),
    ("all4",      "ALL4_highcov_no_REFIND",     8,  "Four\nfocal (8)"),
    ("all4_ref",  "REFIND_plus_ALL4_highcov",   10, "Four focal\n+ donor (10)"),
]

def a_n(n): return sum(1 / i for i in range(1, n))

def read_vcftools(name, n):
    pi = pd.read_csv(STATDIR / f"{name}.windowed.pi", sep="\t")[["BIN_START", "PI"]]
    td = pd.read_csv(STATDIR / f"{name}.Tajima.D", sep="\t")
    td["BIN_START"] = td["BIN_START"] + 1                 # vcftools TajimaD bins are 0-based
    td["theta"] = td["N_SNPS"] / a_n(n) / W
    td["D"] = td["TajimaD"] if n >= 4 else np.nan
    return pi.merge(td[["BIN_START", "theta", "D"]], on="BIN_START", how="inner")

def read_pseudohaploid():
    counts = {}
    with gzip.open(VCFDIR / f"{SAMPLE}_only.snps.biallelic.filtered.vcf.gz", "rt") as fh:
        for line in fh:
            if line[0] == "#": continue
            f = line.split("\t", 10)
            gt = f[9].split(":")[0].replace("|", "/")
            if gt in ("0/1", "1/0", "1/1"):
                b = (int(f[1]) - 1) // W * W + 1
                counts[b] = counts.get(b, 0) + 1
    df = pd.DataFrame({"BIN_START": list(counts), "S": list(counts.values())})
    df["PI"] = df["theta"] = (2 / 3) * df["S"] / W
    df["D"] = np.nan
    return df[["BIN_START", "PI", "theta", "D"]]

data = {}
for key, name, n, _ in CONFIGS:
    data[key] = read_pseudohaploid() if key == "pseudo" else read_vcftools(name, n)
# restrict all configurations to the windows present in every analysis
common = set.intersection(*(set(d["BIN_START"]) for d in data.values()))
for k in data: data[k] = data[k][data[k]["BIN_START"].isin(common)].sort_values("BIN_START")

rows = []
for key, name, n, label in CONFIGS:
    d = data[key]
    rows.append(dict(sample=SAMPLE, srr_id=SRR_ID, configuration=key, haplotypes=n,
                     windows=len(d),
                     pi_mean=d["PI"].mean(), pi_sd=d["PI"].std(),
                     theta_mean=d["theta"].mean(), theta_sd=d["theta"].std(),
                     D_mean=d["D"].mean(), D_sd=d["D"].std()))
summary = pd.DataFrame(rows)
summary.to_csv(PLOTDIR / f"{SRR_ID}_publication_popgen_summary_rev1.csv", index=False)

plt.rcParams.update({"font.family": "sans-serif", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
panels = [("PI", "(a) Nucleotide diversity ($\\pi$)", "$\\pi$ per bp"),
          ("theta", "(b) Watterson's $\\theta_W$", "$\\theta_W$ per bp"),
          ("D", "(c) Tajima's $D$", "Tajima's $D$")]
for ax, (col, title, ylab) in zip(axes, panels):
    vals = [data[k][col].dropna().values for k, *_ in CONFIGS]
    pos = np.arange(len(CONFIGS))
    show = [i for i, v in enumerate(vals) if len(v)]
    bp = ax.boxplot([vals[i] for i in show], positions=pos[show], widths=0.6,
                    patch_artist=True, showfliers=False,
                    medianprops={"color": "black"},
                    boxprops={"facecolor": "#d9d9d9", "edgecolor": "#252525"})
    for i in show:
        ax.plot(pos[i], vals[i].mean(), marker="D", ms=4, color="black", zorder=5)
    for i in set(range(len(CONFIGS))) - set(show):
        ax.text(pos[i], 0.5, "undefined", transform=ax.get_xaxis_transform(),
                ha="center", va="center", fontsize=8, color="#525252", rotation=90)
    if col == "D": ax.axhline(0, color="#969696", lw=0.6, ls=":")
    ax.set_xticks(pos); ax.set_xticklabels([c[3] for c in CONFIGS], fontsize=8, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_title(title, fontsize=10); ax.set_ylabel(ylab)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(PLOTDIR / f"{SRR_ID}_popgen_configurations_rev1.{ext}", dpi=300)
print(summary.to_string(index=False))
