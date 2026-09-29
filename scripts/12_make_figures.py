"""
Renders the 5 figures documenting the cryptic-pocket hypothesis results from
08_correlate_crypticity.py and 11_stratify_by_family.py:

  fig1_ligand_presence_apo_variability.png   hypothesis A (pooled)
  fig2_crypticity_vs_ligand_size.png         hypothesis B (pooled, weak effect)
  fig3_crypticity_vs_diversity.png           hypothesis C (pooled)
  fig4_family_diversity_correlation_top10.png  hypothesis C, per-family "most
                                                evident" (top_n by diversity_rho)
  fig5_crypticity_distribution_by_family_top10.png  family stratification
                                                (top_n by n_sites, most robust)
  fig6_crypticity_vs_affinity.png            hypothesis D (pooled, no effect —
                                                skipped if 13 hasn't been run)

Output: results/figures/*.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from common import load_config

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold",
})
PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
           "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"]


def fmt_p(p: float) -> str:
    """scipy underflows p to exactly 0.0 for very strong/large-n correlations
    (a real float limitation, not a p-value of zero) — report a true upper
    bound instead of the misleading '0.0e+00'."""
    return f"p={p:.1e}" if p > 0 else "p<1e-300"


def savefig(fig, out_dir: Path, name: str):
    path = out_dir / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] {path}")


def fig1_ligand_presence(site_df: pd.DataFrame, out_dir: Path):
    sub = site_df.dropna(subset=["apo_volume_range"])
    with_l = sub[sub["has_ligand"]]["apo_volume_range"]
    without_l = sub[~sub["has_ligand"]]["apo_volume_range"]
    u, p = stats.mannwhitneyu(with_l, without_l, alternative="greater")

    fig, ax = plt.subplots(figsize=(5, 5))
    bp = ax.boxplot([without_l, with_l], labels=["apo-only", "ligand-bound"],
                     showfliers=False, patch_artist=True, widths=0.5)
    for patch, color in zip(bp["boxes"], [PALETTE[7], PALETTE[0]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_yscale("log")
    ax.set_ylabel("apo volume range across crystal forms (A^3, log scale)")
    ax.set_title("A) Apo pocket variability predicts ligand presence")
    ax.text(0.5, 0.02, f"Mann-Whitney U, n={len(with_l)}/{len(without_l)}, p={p:.1e}",
            transform=ax.transAxes, ha="center", fontsize=8, color="dimgray")
    savefig(fig, out_dir, "fig1_ligand_presence_apo_variability.png")


def fig2_ligand_size(crypt: pd.DataFrame, ligands: pd.DataFrame, out_dir: Path):
    def first_code(codes):
        codes = str(codes).split(";") if pd.notna(codes) else []
        return codes[0] if codes else None

    c = crypt.copy()
    c["first_ligand"] = c["holo_ligand_codes_at_max"].apply(first_code)
    c = c.merge(ligands.set_index("het_code")[["heavy_atoms"]],
                left_on="first_ligand", right_index=True, how="left")
    sub = c.dropna(subset=["crypticity_index", "heavy_atoms"])
    rho, p = stats.spearmanr(sub["crypticity_index"], sub["heavy_atoms"])

    fig, ax = plt.subplots(figsize=(5.5, 5))
    hb = ax.hexbin(sub["heavy_atoms"], sub["crypticity_index"], gridsize=40,
                    cmap="Blues", mincnt=1, bins="log")
    fig.colorbar(hb, ax=ax, label="sites (log count)")
    ax.set_xlabel("bound ligand heavy-atom count")
    ax.set_ylabel("crypticity_index")
    ax.set_ylim(0, sub["crypticity_index"].quantile(0.99))
    ax.set_title("B) Ligand size vs. crypticity — weak, if any, effect")
    ax.text(0.97, 0.97, f"Spearman rho={rho:.3f}, p={p:.1e}, n={len(sub)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color="dimgray")
    savefig(fig, out_dir, "fig2_crypticity_vs_ligand_size.png")


def fig3_diversity(crypt: pd.DataFrame, diversity: pd.DataFrame, out_dir: Path):
    sub = diversity.dropna(subset=["crypticity_index", "n_distinct_scaffolds"])
    bins = [0, 1, 2, 5, 10, 30, np.inf]
    labels = ["1", "2", "3-5", "6-10", "11-30", "31+"]
    sub = sub.copy()
    sub["scaffold_bucket"] = pd.cut(sub["n_distinct_scaffolds"], bins=bins, labels=labels)
    rho, p = stats.spearmanr(sub["crypticity_index"], sub["n_distinct_scaffolds"])

    fig, ax = plt.subplots(figsize=(6, 5))
    data = [sub[sub["scaffold_bucket"] == lab]["crypticity_index"].dropna() for lab in labels]
    bp = ax.boxplot(data, labels=labels, showfliers=False, patch_artist=True, widths=0.6)
    for patch in bp["boxes"]:
        patch.set_facecolor(PALETTE[2])
        patch.set_alpha(0.7)
    ax.set_yscale("log")
    ax.set_xlabel("distinct ligand scaffolds observed at this site (bucketed)")
    ax.set_ylabel("crypticity_index (log scale)")
    ax.set_title("C) Chemical diversity vs. crypticity (pooled, all proteins)")
    ax.text(0.02, 0.97, f"Spearman rho={rho:.3f}, {fmt_p(p)}, n={len(sub)}",
            transform=ax.transAxes, ha="left", va="top", fontsize=8, color="dimgray")
    savefig(fig, out_dir, "fig3_crypticity_vs_diversity.png")


def fig6_affinity(affinity_df: pd.DataFrame, out_dir: Path):
    sub = affinity_df.dropna(subset=["crypticity_index", "pchembl_value"])
    rho, p = stats.spearmanr(sub["crypticity_index"], sub["pchembl_value"])

    fig, ax = plt.subplots(figsize=(5.5, 5))
    hb = ax.hexbin(sub["pchembl_value"], sub["crypticity_index"], gridsize=35,
                    cmap="Purples", mincnt=1, bins="log")
    fig.colorbar(hb, ax=ax, label="sites (log count)")
    ax.set_xlabel("ligand potency (PCHEMBL_VALUE, higher = more potent)")
    ax.set_ylabel("crypticity_index")
    ax.set_ylim(0, sub["crypticity_index"].quantile(0.99))
    ax.set_title("D) Binding affinity vs. crypticity — no meaningful effect")
    ax.text(0.97, 0.97, f"Spearman rho={rho:.3f}, {fmt_p(p)}, n={len(sub)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color="dimgray")
    savefig(fig, out_dir, "fig6_crypticity_vs_affinity.png")


def _short_family(name: str, maxlen: int = 38) -> str:
    name = str(name).replace("Belongs to the ", "").replace(" family", "")
    return name if len(name) <= maxlen else name[: maxlen - 1] + "…"


def fig4_family_diversity_top(fam_df: pd.DataFrame, top_n: int, out_dir: Path):
    testable = fam_df.dropna(subset=["diversity_rho"])
    top = testable.sort_values("diversity_rho", ascending=False).head(top_n).iloc[::-1]

    fig, ax = plt.subplots(figsize=(7, 0.5 * len(top) + 1.5))
    y = np.arange(len(top))
    ax.barh(y, top["diversity_rho"], color=PALETTE[3], alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels([_short_family(f) for f in top["protein_family"]], fontsize=8)
    for yi, (rho, p, n) in enumerate(zip(top["diversity_rho"], top["diversity_p"], top["n_sites"])):
        ax.text(rho + 0.01, yi, f"p={p:.0e}, n={n}", va="center", fontsize=7, color="dimgray")
    ax.set_xlabel("within-family Spearman rho (crypticity vs. n distinct scaffolds)")
    ax.set_title(f"C, stratified) Top {top_n} families — strongest within-family signal")
    ax.set_xlim(0, min(1.0, top["diversity_rho"].max() * 1.35))
    savefig(fig, out_dir, "fig4_family_diversity_correlation_top10.png")


def fig5_family_distribution_top(crypt: pd.DataFrame, fam_df: pd.DataFrame, top_n: int, out_dir: Path):
    top = fam_df.sort_values("n_sites", ascending=False).head(top_n)
    families = top["protein_family"].tolist()
    data = [crypt[crypt["protein_family"] == f]["crypticity_index"].values for f in families]
    groups_for_test = [crypt[crypt["protein_family"] == f]["crypticity_index"].values
                        for f in fam_df["protein_family"]]
    h, p = stats.kruskal(*groups_for_test)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    bp = ax.boxplot(data, labels=[_short_family(f, 24) for f in families],
                     showfliers=False, patch_artist=True, widths=0.6)
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(PALETTE[i % len(PALETTE)])
        patch.set_alpha(0.75)
    ax.set_yscale("log")
    ax.set_ylabel("crypticity_index (log scale)")
    ax.set_title(f"Crypticity varies by protein family — top {top_n} by n_sites")
    ax.text(0.99, 0.97, f"Kruskal-Wallis across all {len(fam_df)} families\nwith "
            f">= {int(fam_df['n_sites'].min())} sites: H={h:.0f}, {fmt_p(p)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color="dimgray")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", fontsize=8)
    savefig(fig, out_dir, "fig5_crypticity_distribution_by_family_top10.png")


def main():
    cfg = load_config()
    strat = cfg["stratification"]
    out_dir = Path(cfg["paths"]["results"]) / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        crypt = pd.read_csv(f"{cfg['paths']['results']}/crypticity_index_with_family.csv")
        diversity = pd.read_csv(f"{cfg['paths']['results']}/diversity_analysis.csv")
        site_df = pd.read_csv(f"{cfg['paths']['results']}/site_level_analysis.csv")
        ligands = pd.read_csv(f"{cfg['paths']['ligands']}/ligand_descriptors.csv")
        fam_df = pd.read_csv(f"{cfg['paths']['results']}/family_stratification.csv")
    except FileNotFoundError as e:
        sys.exit(f"Missing input ({e}). Run 06/07/08/10/11 first.")

    try:
        affinity_df = pd.read_csv(f"{cfg['paths']['results']}/crypticity_index_with_affinity.csv")
    except FileNotFoundError:
        affinity_df = None
        print("[skip] fig6: missing results/crypticity_index_with_affinity.csv. Run 13 first.")

    diversity = diversity.merge(
        crypt[["pocket_site_id", "protein_family"]], on="pocket_site_id", how="left")

    fig1_ligand_presence(site_df, out_dir)
    fig2_ligand_size(crypt, ligands, out_dir)
    fig3_diversity(crypt, diversity, out_dir)
    fig4_family_diversity_top(fam_df, strat["top_n_families"], out_dir)
    fig5_family_distribution_top(crypt, fam_df, strat["top_n_families"], out_dir)
    n_figs = 5
    if affinity_df is not None:
        fig6_affinity(affinity_df, out_dir)
        n_figs = 6

    print(f"\nWrote {n_figs} figures -> {out_dir}")


if __name__ == "__main__":
    main()
