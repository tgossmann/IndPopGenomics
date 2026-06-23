#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    sys.exit("Usage: python3 plot_triallelic_barplot_highcov.py <OUTDIR>")

OUTDIR = Path(sys.argv[1])
TRIDIR = OUTDIR / "triallelic"
PLOTDIR = OUTDIR / "plots_publication"
PLOTDIR.mkdir(exist_ok=True)

counts_file = TRIDIR / "strict_triallelic_counts.tsv"

df = pd.read_csv(counts_file, sep="\t")

label_map = {
    "H1_only": "Focal 1",
    "H2_only": "Focal 2",
    "H3_only": "Focal 3",
    "H4_only": "Focal 4",
    "REFIND_only": "Reference",
    "REFIND_plus_H1": "Reference + focal 1",
    "REFIND_plus_H2": "Reference + focal 2",
    "REFIND_plus_H3": "Reference + focal 3",
    "REFIND_plus_H4": "Reference + focal 4",
    "ALL4_highcov_no_REFIND": "Four individuals",
    "REFIND_plus_ALL4_highcov": "Four individuals + reference",
}

order = [
    "H1_only",
    "REFIND_plus_H1",
    "H2_only",
    "REFIND_plus_H2",
    "H3_only",
    "REFIND_plus_H3",
    "H4_only",
    "REFIND_plus_H4",
    "REFIND_only",
    "ALL4_highcov_no_REFIND",
    "REFIND_plus_ALL4_highcov",
]

df = df[df["analysis"].isin(order)].copy()
df["label"] = df["analysis"].map(label_map)
df["order"] = df["analysis"].apply(order.index)
df = df.sort_values("order")

plt.figure(figsize=(10, 5.5))
bars = plt.bar(df["label"], df["strict_triallelic_snps"])

plt.ylabel("Strict triallelic SNPs")
plt.xlabel("")
plt.title("Reference donor inclusion reveals additional triallelic SNPs")
plt.xticks(rotation=45, ha="right")

for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        height,
        f"{int(height)}",
        ha="center",
        va="bottom",
        fontsize=8,
    )

plt.tight_layout()

outfile = PLOTDIR / "strict_triallelic_snps_publication_barplot.png"
plt.savefig(outfile, dpi=300)
plt.close()

summary_out = PLOTDIR / "strict_triallelic_snps_publication_barplot_data.csv"
df.to_csv(summary_out, index=False)

print(f"Wrote: {outfile}")
print(f"Wrote: {summary_out}")
