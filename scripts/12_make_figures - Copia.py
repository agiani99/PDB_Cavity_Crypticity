from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.ticker import FixedLocator, LogFormatterMathtext, NullFormatter

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import load_config

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold"
})

PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
           "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"]


def savefig(fig, out_dir: Path, name: str):
    path = out_dir / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] {path}")


def family_label(name: str) -> str:
    """Preserve complete family names; insert line breaks only after periods."""
    name = str(name).replace("Belongs to the ", "").replace(" family", "")
    return " ".join(name.split()).replace(".", ".\n")


def fig1_ligand_presence(site_df: pd.DataFrame, out_dir: Path):
    sub = site_df.dropna(subset=["apo_volume_range"])
    with_l = sub[sub["has_ligand"]]["apo_volume_range"]
    without_l = sub[~sub["has_ligand"]]["apo_volume_range"]

    fig, ax = plt.subplots(figsize=(5, 5))
    bp = ax.boxplot([without_l, with_l], tick_labels=["Apo-only", "Ligand-bound"],
                    showfliers=False, patch_artist=True, widths=0.5)

    for patch, color in zip(bp["boxes"], [PALETTE[7], PALETTE[0]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_yscale("log")
    ax.set_ylabel(r"Apo pocket volume range ($\mathrm{\AA}^3$)")
    savefig(fig, out_dir, "fig2_ligand_presence_apo_variability.png")


def fig2_ligand_size(crypt: pd.DataFrame, ligands: pd.DataFrame, out_dir: Path):
    def first_code(codes):
        codes = str(codes).split(";") if pd.notna(codes) else []
        return codes[0] if codes else None

    c = crypt.copy()
    c["first_ligand"] = c["holo_ligand_codes_at_max"].apply(first_code)
    c = c.merge(ligands.set_index("het_code")[["heavy_atoms"]],
                left_on="first_ligand", right_index=True, how="left")
    sub = c.dropna(subset=["crypticity_index", "heavy_atoms"])

    fig, ax = plt.subplots(figsize=(5.5, 5))
    hb = ax.hexbin(sub["heavy_atoms"], sub["crypticity_index"],
                   gridsize=40, cmap="Blues", mincnt=1, bins="log")
    fig.colorbar(hb, ax=ax, label="sites (log count)")
    ax.set_xlabel("bound ligand heavy-atom count")
    ax.set_ylabel("crypticity index")
    ax.set_ylim(0, sub["crypticity_index"].quantile(0.99))
    savefig(fig, out_dir, "fig3_crypticity_vs_ligand_size.png")


def fig3_diversity(diversity: pd.DataFrame, out_dir: Path):
    sub = diversity.dropna(subset=["crypticity_index", "n_distinct_scaffolds"])
    bins = [0, 1, 2, 5, 10, 30, np.inf]
    labels = ["1", "2", "3-5", "6-10", "11-30", "31+"]

    sub = sub.copy()
    sub["scaffold_bucket"] = pd.cut(sub["n_distinct_scaffolds"], bins=bins, labels=labels)

    fig, ax = plt.subplots(figsize=(6, 5))
    data = [sub[sub["scaffold_bucket"] == lab]["crypticity_index"].dropna() for lab in labels]
    bp = ax.boxplot(data, tick_labels=labels, showfliers=False,
                    patch_artist=True, widths=0.6)

    for patch in bp["boxes"]:
        patch.set_facecolor(PALETTE[2])
        patch.set_alpha(0.7)

    ax.set_yscale("log")
    ax.set_xlabel("distinct ligand scaffolds observed at this site (bucketed)")
    ax.set_ylabel("crypticity index (log scale)")
    savefig(fig, out_dir, "fig4_crypticity_vs_diversity.png")


def fig4_family_diversity_top(fam_df: pd.DataFrame, top_n: int, out_dir: Path):
    testable = fam_df.dropna(subset=["diversity_rho"])
    top = testable.sort_values("diversity_rho", ascending=False).head(top_n).iloc[::-1]

    fig, ax = plt.subplots(figsize=(7, 0.5 * len(top) + 1.5))
    y = np.arange(len(top))
    ax.barh(y, top["diversity_rho"], color=PALETTE[3], alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels([family_label(f) for f in top["protein_family"]], fontsize=8)

    for yi, (rho, p, n) in enumerate(zip(top["diversity_rho"], top["diversity_p"], top["n_sites"])):
        ax.text(rho + 0.01, yi, f"n={int(n)}", va="center", ha="left",
                fontsize=7, color="dimgray")

    ax.set_xlabel("within-family Spearman rho (crypticity vs. n distinct scaffolds)")
    max_rho = top["diversity_rho"].max()
    ax.set_xlim(0, min(1.0, max_rho * 1.35))
    fig.subplots_adjust(left=0.35)
    savefig(fig, out_dir, "fig5_family_diversity_correlation_top10.png")


def fig5_family_distribution_top(crypt: pd.DataFrame, fam_df: pd.DataFrame,
                                 top_n: int, out_dir: Path):
    top = fam_df.sort_values("n_sites", ascending=False).head(top_n)
    families = top["protein_family"].tolist()
    data = [crypt[crypt["protein_family"] == f]["crypticity_index"].values for f in families]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    bp = ax.boxplot(data, tick_labels=[family_label(f) for f in families],
                    showfliers=False, patch_artist=True, widths=0.6)

    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(PALETTE[i % len(PALETTE)])
        patch.set_alpha(0.75)

    ax.set_yscale("log")
    ax.yaxis.set_major_locator(FixedLocator([1]))
    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylabel("crypticity index (log scale)")

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor", fontsize=8)
    fig.subplots_adjust(bottom=0.38)
    savefig(fig, out_dir, "fig6_crypticity_distribution_by_family_top10.png")


def fig6_affinity(affinity_df: pd.DataFrame, out_dir: Path):
    sub = affinity_df.dropna(subset=["crypticity_index", "pchembl_value"])

    fig, ax = plt.subplots(figsize=(5.5, 5))
    hb = ax.hexbin(sub["pchembl_value"], sub["crypticity_index"],
                   gridsize=35, cmap="Purples", mincnt=1, bins="log")
    fig.colorbar(hb, ax=ax, label="sites (log count)")
    ax.set_xlabel("ligand potency (pChEMBL)")
    ax.set_ylabel("crypticity index")
    ax.set_ylim(0, sub["crypticity_index"].quantile(0.99))
    savefig(fig, out_dir, "fig7_crypticity_vs_affinity.png")


def fig7_apo_druggability(drug_df: pd.DataFrame, out_dir: Path):
    cryptic = drug_df[drug_df["crypticity_index"] >= 5]["druggability_apo_min_v"].dropna()
    noncryptic = drug_df[drug_df["crypticity_index"] < 2]["druggability_apo_min_v"].dropna()

    fig, ax = plt.subplots(figsize=(6, 5))

    for data, label, color in [
        (noncryptic, f"minimally cryptic (CI<2), n={len(noncryptic)}", PALETTE[0]),
        (cryptic, f"highly cryptic (CI≥5), n={len(cryptic)}", PALETTE[3])
    ]:
        x = np.sort(data.values)
        if len(x) == 0:
            continue
        y = np.arange(1, len(x) + 1) / len(x)
        ax.plot(x, y, label=label, color=color, linewidth=2)

    ax.set_xlabel("apo-state fpocket druggability score")
    ax.set_ylabel("cumulative fraction of sites")
    ax.set_xlim(0, 1)
    ax.legend(loc="lower right", fontsize=8)
    savefig(fig, out_dir, "fig8_apo_druggability_vs_crypticity.png")


def load_csv(path: Path, figure_name: str):
    """Load a CSV independently; return None if it is missing or empty."""
    if not path.exists():
        print(f"[skip] {figure_name}: missing input file: {path}")
        return None
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        print(f"[skip] {figure_name}: input file is empty: {path}")
        return None


def main():
    cfg = load_config()
    strat = cfg["stratification"]
    results_dir = Path(cfg["paths"]["results"])
    ligands_dir = Path(cfg["paths"]["ligands"])
    out_dir = results_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    crypt = load_csv(results_dir / "crypticity_index_with_family.csv", "crypticity data")
    diversity = load_csv(results_dir / "diversity_analysis.csv", "fig3")
    site_df = load_csv(results_dir / "site_level_analysis.csv", "fig1")
    ligands = load_csv(ligands_dir / "ligand_descriptors.csv", "fig2")
    fam_df = load_csv(results_dir / "family_stratification.csv", "fig4 and fig5")
    affinity_df = load_csv(results_dir / "crypticity_index_with_affinity.csv", "fig6")
    drug_df = load_csv(results_dir / "druggability_analysis.csv", "fig7")

    n_figs = 0

    if site_df is not None:
        fig1_ligand_presence(site_df, out_dir)
        n_figs += 1

    if crypt is not None and ligands is not None:
        fig2_ligand_size(crypt, ligands, out_dir)
        n_figs += 1
    else:
        print("[skip] fig2: requires crypticity_index_with_family.csv and ligand_descriptors.csv.")

    if diversity is not None:
        fig3_diversity(diversity, out_dir)
        n_figs += 1

    if fam_df is not None:
        fig4_family_diversity_top(fam_df, strat["top_n_families"], out_dir)
        n_figs += 1

    if crypt is not None and fam_df is not None:
        fig5_family_distribution_top(crypt, fam_df, strat["top_n_families"], out_dir)
        n_figs += 1
    else:
        print("[skip] fig5: requires crypticity_index_with_family.csv and family_stratification.csv.")

    if affinity_df is not None:
        fig6_affinity(affinity_df, out_dir)
        n_figs += 1

    if drug_df is not None:
        fig7_apo_druggability(drug_df, out_dir)
        n_figs += 1

    print(f"\nWrote {n_figs} figures -> {out_dir}")


if __name__ == "__main__":
    main()
