#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


if len(sys.argv) != 4:
    sys.exit(
        "Usage: python3 plot_publication_popgen_highcov.py "
        "<OUTDIR> <H1|H2|H3|H4> <SRR_ID>"
    )

OUTDIR = Path(sys.argv[1])
SAMPLE = sys.argv[2]
SRR_ID = sys.argv[3]

STATDIR = OUTDIR / "stats"
PLOTDIR = OUTDIR / "plots_publication"
PLOTDIR.mkdir(exist_ok=True)

WINDOW_SIZE = 50_000


ANALYSES = {
    "focal_ref": {
        "name": f"REFIND_plus_{SAMPLE}",
        "label": "Focal + Reference (4 haplotypes)",
    },
    "all4": {
        "name": "ALL4_highcov_no_REFIND",
        "label": "Four individuals (8 haplotypes)",
    },
    "all4_ref": {
        "name": "REFIND_plus_ALL4_highcov",
        "label": "Four individuals + Reference (10 haplotypes)",
    },
}


def harmonic_number(n):
    return sum(1 / i for i in range(1, n))


def n_chromosomes(analysis_name):
    if analysis_name.startswith("REFIND_plus_H"):
        return 4
    if analysis_name == "ALL4_highcov_no_REFIND":
        return 8
    if analysis_name == "REFIND_plus_ALL4_highcov":
        return 10
    raise ValueError(f"Unknown analysis: {analysis_name}")


def read_metric(analysis_name, metric):
    if metric == "pi":
        path = STATDIR / f"{analysis_name}.windowed.pi"
        df = pd.read_csv(path, sep="\t").sort_values("BIN_START")
        df["value"] = df["PI"]
        return df

    path = STATDIR / f"{analysis_name}.Tajima.D"
    df = pd.read_csv(path, sep="\t").sort_values("BIN_START")

    if metric == "theta":
        n = n_chromosomes(analysis_name)
        a_n = harmonic_number(n)
        df["value"] = df["N_SNPS"] / a_n / WINDOW_SIZE
        return df

    if metric == "tajimaD":
        df["value"] = df["TajimaD"]
        return df

    raise ValueError(metric)


def add_mean_line(ax, y, color, label):
    ax.axhline(
        y,
        linestyle="--",
        linewidth=1.6,
        alpha=0.85,
        color=color,
        label=label,
    )


def metric_label(metric):
    return {
        "pi": "Nucleotide diversity (π)",
        "theta": "Watterson's θ per bp",
        "tajimaD": "Tajima's D",
    }[metric]


def metric_filename(metric):
    return {
        "pi": "pi",
        "theta": "thetaW",
        "tajimaD": "tajimaD",
    }[metric]


def plot_metric(metric):
    plot_keys = ["focal_ref", "all4", "all4_ref"]

    data = {
        key: read_metric(ANALYSES[key]["name"], metric)
        for key in plot_keys
    }

    means = {
        key: df["value"].mean()
        for key, df in data.items()
    }

    fig, ax = plt.subplots(figsize=(11, 4.8))

    for key in plot_keys:
        info = ANALYSES[key]
        df = data[key]

        line = ax.plot(
            df["BIN_START"],
            df["value"],
            label=info["label"],
        )[0]

        add_mean_line(
            ax,
            means[key],
            line.get_color(),
            f"{info['label']} mean = {means[key]:.6g}",
        )

    ax.set_xlabel("Genomic position")
    ax.set_ylabel(metric_label(metric))
    ax.set_title(metric_label(metric))
    ax.legend(fontsize=8, frameon=False)

    plt.tight_layout()

    outfile = (
        PLOTDIR
        / f"{SRR_ID}_{metric_filename(metric)}_publication_comparison_no_single_individual.png"
    )
    plt.savefig(outfile, dpi=300)
    plt.close()

    row = {
        "sample": SAMPLE,
        "srr_id": SRR_ID,
        "metric": metric,
    }

    for key in plot_keys:
        row[f"{key}_mean"] = means[key]

    return row


rows = []

for metric in ["pi", "theta", "tajimaD"]:
    print(f"Plotting {metric}")
    rows.append(plot_metric(metric))

summary = pd.DataFrame(rows)
summary_file = (
    PLOTDIR
    / f"{SRR_ID}_publication_popgen_summary_no_single_individual.csv"
)
summary.to_csv(summary_file, index=False)

print(f"Done. Output written to: {PLOTDIR.resolve()}")
print(summary)
