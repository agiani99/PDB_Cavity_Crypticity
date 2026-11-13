"""
Alternative Figure 1 as violin plots (replaces the box plot of 12_make_figures.py fig1; does not modify it).

Panel A  distribution of the crypticity index over all scored sites (log10 scale violin, median and mean marked,
         CI = 1 dashed). This is the distribution the Results text describes ("right-skewed, median 1.542, mean 1.962").
Panel B  Hypothesis A: apo volume range across structures for sites that later bind a ligand vs apo-only sites
         (same data and Mann-Whitney test as the old Figure 1).

Inputs : results/crypticity_index.csv            (panel A)
         results/site_level_analysis.csv         (panel B)
Output : results/figures/fig1_violin_distribution_and_hypA.png  (300 dpi)
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import load_config

BLUE, GREY, ORANGE = "#5F7392", "#9AA3AD", "#E08A2E"


def violin(ax, data, positions, colors, log10=True):
    """Violin on log10 values with quartile box; returns the plotted arrays."""
    vals = [np.log10(np.asarray(d, float)) if log10 else np.asarray(d, float) for d in data]
    parts = ax.violinplot(vals, positions=positions, widths=0.8, showextrema=False, showmedians=False)
    for body, col in zip(parts["bodies"], colors):
        body.set_facecolor(col)
        body.set_edgecolor("#33425B")
        body.set_alpha(0.75)
    for pos, v in zip(positions, vals):
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.vlines(pos, q1, q3, color="#33425B", lw=4, zorder=3)
        ax.scatter([pos], [med], color="white", s=14, zorder=4)
    return vals


def main():
    cfg = load_config()
    res = Path(cfg["paths"]["results"])
    ci = pd.read_csv(res / "crypticity_index.csv")["crypticity_index"].dropna()
    ci = ci[ci > 0]
    site = pd.read_csv(res / "site_level_analysis.csv").dropna(subset=["apo_volume_range"])
    site = site[site["apo_volume_range"] > 0]
    with_l, without_l = site[site.has_ligand].apo_volume_range, site[~site.has_ligand].apo_volume_range
    p = stats.mannwhitneyu(with_l, without_l, alternative="greater").pvalue

    fig, (a, b) = plt.subplots(1, 2, figsize=(8.2, 4.4), gridspec_kw={"width_ratios": [1, 1.25]})

    violin(a, [ci], [1], [BLUE])
    a.axhline(0, color="k", lw=0.8, ls="--")
    a.scatter([1], [np.log10(ci.mean())], marker="D", color=ORANGE, s=26, zorder=5, label=f"mean {ci.mean():.2f}")
    a.scatter([1], [np.log10(ci.median())], marker="o", color="white", edgecolor="#33425B", s=30, zorder=5,
              label=f"median {ci.median():.2f}")
    ticks = [0.1, 0.3, 1, 3, 10, 30, 100]
    a.set_yticks(np.log10(ticks))
    a.set_yticklabels([str(t) for t in ticks])
    a.set_xticks([1]); a.set_xticklabels([f"all sites\n(n={len(ci):,})"])
    a.set_xlim(0.4, 1.6)
    a.set_ylabel("crypticity index, CI (log scale)")
    a.set_title("A  Distribution of CI", loc="left", fontsize=10, fontweight="bold")
    a.text(1.55, 0.04, "CI = 1", fontsize=7, ha="right", va="bottom", color="dimgray")
    a.legend(loc="upper right", fontsize=8, frameon=False)
    a.text(0.03, 0.03, f"CI < 1: {100 * (ci < 1).mean():.1f}%\nCI > 2: {100 * (ci > 2).mean():.1f}%\nCI > 5: {100 * (ci > 5).mean():.1f}%",
           transform=a.transAxes, fontsize=7.5, va="bottom", ha="left", color="#33425B")

    violin(b, [without_l, with_l], [1, 2], [GREY, BLUE])
    ticks = [1, 10, 100, 1000, 10000]
    b.set_yticks(np.log10(ticks))
    b.set_yticklabels([f"{t:,}" for t in ticks])
    b.set_ylim(top=np.log10(60000))   # headroom so the annotation does not overlap the violin tips
    b.set_xticks([1, 2])
    b.set_xticklabels([f"apo-only sites\n(n={len(without_l):,})", f"ligand-bound sites\n(n={len(with_l):,})"])
    b.set_ylabel("apo volume range across structures (Å³, log scale)")
    b.set_title("B  Apo pocket variability (Hypothesis A)", loc="left", fontsize=10, fontweight="bold")
    b.text(0.5, 0.97, f"medians {without_l.median():.0f} vs {with_l.median():.0f} Å³; Mann–Whitney p = {p:.1e}",
           transform=b.transAxes, ha="center", va="top", fontsize=8, color="#33425B")

    for ax in (a, b):
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = res / "figures" / "fig1_violin_distribution_and_hypA.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    print(f"wrote {out}\nCI: n={len(ci)}, median={ci.median():.3f}, mean={ci.mean():.3f}; "
          f"apo range medians {without_l.median():.1f} (n={len(without_l)}) vs {with_l.median():.1f} (n={len(with_l)}), p={p:.2e}")


if __name__ == "__main__":
    main()
