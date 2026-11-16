
"""
Alternative Figure 1: violin plot of the crypticity-index distribution.

Panel A only: distribution of the crypticity index over all scored sites,
shown as a violin plot on a log10 scale, with median and mean marked
and CI = 1 indicated by a dashed reference line.

This is an alternative to Panel A of Figure 1 in 12_make_figures.py
and does not modify that script.

Input:
    results/crypticity_index.csv

Output:
    results/figures/fig1_violin_distribution.png (300 dpi)
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import load_config

BLUE, ORANGE = "#5F7392", "#E08A2E"


def violin(ax, data, position=1, color=BLUE, log10=True):
    """Draw a violin on log10 values with quartile box; return plotted values."""
    values = np.asarray(data, float)
    vals = np.log10(values) if log10 else values

    parts = ax.violinplot(
        [vals],
        positions=[position],
        widths=0.8,
        showextrema=False,
        showmedians=False
    )

    for body in parts["bodies"]:
        body.set_facecolor(color)
        body.set_edgecolor("#33425B")
        body.set_alpha(0.75)

    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    ax.vlines(position, q1, q3, color="#33425B", lw=4, zorder=3)
    ax.scatter([position], [med], color="white", s=14, zorder=4)

    return vals


def main():
    cfg = load_config()
    res = Path(cfg["paths"]["results"])

    ci = pd.read_csv(res / "crypticity_index.csv")["crypticity_index"].dropna()
    ci = ci[ci > 0]

    fig, ax = plt.subplots(1, 1, figsize=(4.2, 4.4))

    violin(ax, ci, position=1, color=BLUE)

    ax.axhline(0, color="k", lw=0.8, ls="--")

    ax.scatter(
        [1],
        [np.log10(ci.mean())],
        marker="D",
        color=ORANGE,
        s=26,
        zorder=5,
        label=f"mean {ci.mean():.2f}"
    )

    ax.scatter(
        [1],
        [np.log10(ci.median())],
        marker="o",
        color="white",
        edgecolor="#33425B",
        s=30,
        zorder=5,
        label=f"median {ci.median():.2f}"
    )

    ticks = [0.1, 0.3, 1, 3, 10, 30, 100]
    ax.set_yticks(np.log10(ticks))
    ax.set_yticklabels([str(t) for t in ticks])

    ax.set_xticks([1])
    ax.set_xticklabels([f"all sites (n={len(ci):,})"])
    ax.set_xlim(0.4, 1.6)

    ax.set_ylabel("crypticity index (log scale)")

    ax.text(
        1.55,
        0.04,
        "CI = 1",
        fontsize=7,
        ha="right",
        va="bottom",
        color="dimgray"
    )


    ax.text(
        0.03,
        0.03,
        f"CI < 1: {100 * (ci < 1).mean():.1f}%\n"
        f"CI > 2: {100 * (ci > 2).mean():.1f}%\n"
        f"CI > 5: {100 * (ci > 5).mean():.1f}%",
        transform=ax.transAxes,
        fontsize=7.5,
        va="bottom",
        ha="left",
        color="#33425B"
    )

    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()

    out = res / "figures" / "fig1_violin_distribution.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")

    print(
        f"wrote {out}\n"
        f"CI: n={len(ci)}, median={ci.median():.3f}, "
        f"mean={ci.mean():.3f}"
    )

    plt.close(fig)


if __name__ == "__main__":
    main()
