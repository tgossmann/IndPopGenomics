#!/usr/bin/env python3
"""Figure 2: observable spectra with the reference donor resolved vs. not resolved.
Reads sim_summary.csv (from summarise_sim.py); writes SFS_mosaic_simulation.pdf."""
import csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})
rows = {r["scenario"]: r for r in csv.DictReader(open("sim_summary.csv"))}
L, MU, NE = 2_000_000, 1e-8, 10_000
theta = 4 * NE * MU * L                      # expected theta per replicate (constant size)
expect = {"M1": 4/3*theta, "M2": theta/2, "K1": theta, "K2": theta/2, "lost": theta/3}
titles = {"constant": "(a) Constant size", "expansion": "(b) Expansion", "contraction": "(c) Contraction"}
keys   = ["M1", "M2", "K1", "K2", "lost"]
labels = ["1 of 4", "2 of 4", "1 of 3", "2 of 3", "lost"]
xpos   = np.array([0, 1, 2.6, 3.6, 4.6])
fill   = ["#bdbdbd", "#525252", "#bdbdbd", "#525252", "white"]
fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
for ax, sc in zip(axes, titles):
    r = rows[sc]
    D = float(r["D4_mean"]); lostpct = 100 * float(r["lost_frac"])
    for i, k in enumerate(keys):
        ax.bar(xpos[i], float(r[k + "_mean"]), width=0.8, color=fill[i], edgecolor="#252525",
               hatch="///" if k == "lost" else None, linewidth=0.8,
               yerr=float(r[k + "_sd"]), capsize=3, error_kw={"lw": 0.8})
        if sc == "constant":
            ax.plot(xpos[i], expect[k], marker="D", ms=5, color="black", zorder=5)
    ax.set_xticks(xpos); ax.set_xticklabels(labels)
    ax.axvline(1.8, color="#969696", lw=0.6, ls=":")
    ymax = 1800
    ax.text(0.5, ymax * 0.97, "donor resolved\n(4 haplotypes)", ha="center", va="top", fontsize=8.5)
    ax.text(3.6, ymax * 0.97, "donor not resolved\n(3 sequences)", ha="center", va="top", fontsize=8.5)
    ax.set_title(f"{titles[sc]}\nTajima's D (4 haplotypes) = {D:.2f}; lost {lostpct:.0f}%", fontsize=9.5)
    ax.set_ylim(0, ymax)
axes[0].set_ylabel("Segregating sites per 2 Mb")
fig.tight_layout()
fig.savefig("SFS_mosaic_simulation.pdf"); fig.savefig("SFS_mosaic_simulation.png", dpi=130)
