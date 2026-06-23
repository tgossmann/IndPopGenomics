#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


if len(sys.argv) != 2:
    sys.exit("Usage: python3 plot_highcov_popgen_comparison.py <OUTDIR>")

OUTDIR = Path(sys.argv[1])
STATDIR = OUTDIR / "stats"
PLOTDIR = OUTDIR / "plots_popgen_comparison"
PLOTDIR.mkdir(exist_ok=True)

WINDOW_SIZE = 50_000


def harmonic_number(n):
    return sum(1 / i for i in range(1, n))


def n_chromosomes(name):
    if name.endswith("_only"):
        return 2
    if name.startswith("REFIND_plus_H"):
        return 4
    if name == "ALL4_highcov_no_REFIND":
        return 8
    if name == "REFIND_plus_ALL4_highcov":
        return 10
    raise ValueError(f"Unknown analysis: {name}")


def read_pi(name):
    path = STATDIR / f"{name}.windowed.pi"
    df = pd.read_csv(path, sep="\t")
    df = df.sort_values("BIN_START")
    return df


def read_taj_theta(name):
    path = STATDIR / f"{name}.Tajima.D"
    df = pd.read_csv(path, sep="\t")
    df = df.sort_values("BIN_START")

    n = n_chromosomes(name)
    a_n = harmonic_number(n)

    df["theta_W"] = df["N_SNPS"] / a_n / WINDOW_SIZE
    df["n_chromosomes"] = n
    return df


def add_mean(ax, y, color, label):
    ax.axhline(
        y,
        linestyle="--",
        linewidth=1.6,
        alpha=0.85,
        color=color,
        label=label,
    )


def plot_metric_for_sample(sample, metric):
    refplus_name = f"REFIND_plus_{sample}"
    all4_name = "ALL4_highcov_no_REFIND"
    all4ref_name = "REFIND_plus_ALL4_highcov"

    if metric == "pi":
        refplus = read_pi(refplus_name)
        all4 = read_pi(all4_name)
        all4ref = read_pi(all4ref_name)
        value_col = "PI"
        ylabel = "Windowed π"
        title_metric = "π"
    elif metric == "theta":
        refplus = read_taj_theta(refplus_name)
        all4 = read_taj_theta(all4_name)
        all4ref = read_taj_theta(all4ref_name)
        value_col = "theta_W"
        ylabel = "Watterson's θ per bp"
        title_metric = "θW"
    elif metric == "tajima":
        refplus = read_taj_theta(refplus_name)
        all4 = read_taj_theta(all4_name)
        all4ref = read_taj_theta(all4ref_name)
        value_col = "TajimaD"
        ylabel = "Tajima's D"
        title_metric = "Tajima's D"
    else:
        raise ValueError(metric)

    merged = refplus.merge(
        all4,
        on=["CHROM", "BIN_START"],
        suffixes=("_refplus", "_all4"),
    )

    merged[f"delta_{metric}"] = merged[f"{value_col}_refplus"] - merged[f"{value_col}_all4"]

    mean_refplus = refplus[value_col].mean()
    mean_all4 = all4[value_col].mean()
    mean_all4ref = all4ref[value_col].mean()
    mean_delta = merged[f"delta_{metric}"].mean()

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

    l1 = axes[0].plot(
        refplus["BIN_START"],
        refplus[value_col],
        label=f"REFIND + {sample} (4 haplotypes)",
    )[0]

    l2 = axes[0].plot(
        all4["BIN_START"],
        all4[value_col],
        label="ALL4 no REFIND (8 haplotypes)",
    )[0]

    l3 = axes[0].plot(
        all4ref["BIN_START"],
        all4ref[value_col],
        label="REFIND + ALL4 (10 haplotypes)",
    )[0]

    add_mean(axes[0], mean_refplus, l1.get_color(), f"REFIND+{sample} mean = {mean_refplus:.6g}")
    add_mean(axes[0], mean_all4, l2.get_color(), f"ALL4 mean = {mean_all4:.6g}")
    add_mean(axes[0], mean_all4ref, l3.get_color(), f"REFIND+ALL4 mean = {mean_all4ref:.6g}")

    axes[0].set_ylabel(ylabel)
    axes[0].set_title(f"{title_metric}: REFIND + {sample} versus population panels")
    axes[0].legend(fontsize=8)

    d = axes[1].plot(
        merged["BIN_START"],
        merged[f"delta_{metric}"],
        label=f"REFIND+{sample} minus ALL4",
    )[0]

    axes[1].axhline(0, linestyle="--", color="black", linewidth=1, alpha=0.7)
    add_mean(axes[1], mean_delta, d.get_color(), f"mean difference = {mean_delta:.6g}")

    axes[1].set_xlabel("Genomic position")
    axes[1].set_ylabel(f"Δ {title_metric}")
    axes[1].set_title(f"Difference: REFIND + {sample} minus ALL4 no REFIND")
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    outfile = PLOTDIR / f"{sample}_{metric}_refind_plus_sample_vs_all4.png"
    plt.savefig(outfile, dpi=300)
    plt.close()

    return {
        "sample": sample,
        "metric": metric,
        "mean_REFIND_plus_sample": mean_refplus,
        "mean_ALL4_no_REFIND": mean_all4,
        "mean_REFIND_plus_ALL4": mean_all4ref,
        "mean_REFIND_plus_sample_minus_ALL4": mean_delta,
    }


rows = []

for sample in ["H1", "H2", "H3", "H4"]:
    for metric in ["pi", "theta", "tajima"]:
        print(f"Plotting {sample} {metric}")
        rows.append(plot_metric_for_sample(sample, metric))

summary = pd.DataFrame(rows)
summary.to_csv(PLOTDIR / "popgen_summary_refind_sample_vs_all4.csv", index=False)

print(f"Done. Plots written to: {PLOTDIR.resolve()}")
print(summary)
