"""
Formalizes the ad hoc family-stratification check: does crypticity (and the
crypticity-vs-diversity relationship from 08) hold up within individual
UniProt protein families, or is the pooled result just a between-family
artifact (e.g. one oversampled family dominating the pooled numbers)?

Only families with >= config.yaml [stratification].min_family_n scored
sites are tested individually — smaller families still count toward the
overall Kruskal-Wallis but don't get their own per-family row (too few
sites to say anything meaningful).

Output:
    results/family_stratification.csv
        protein_family, n_sites, median_crypticity, mean_crypticity,
        diversity_rho, diversity_p
    results/family_stratification.txt
        overall Kruskal-Wallis result + the top_n_families "most evident"
        by diversity_rho and by n_sites (used to pick panels for figures)
"""
import sys
from pathlib import Path

import pandas as pd
from scipy import stats
from common import load_config


def main():
    cfg = load_config()
    strat = cfg["stratification"]
    try:
        crypt = pd.read_csv(f"{cfg['paths']['results']}/crypticity_index_with_family.csv")
        diversity = pd.read_csv(f"{cfg['paths']['results']}/diversity_analysis.csv")
    except FileNotFoundError as e:
        sys.exit(f"Missing input ({e}). Run 06/08/10 first.")

    diversity = diversity.merge(
        crypt[["pocket_site_id", "protein_family"]], on="pocket_site_id", how="left")

    fam_counts = crypt.groupby("protein_family").size()
    big_families = fam_counts[fam_counts >= strat["min_family_n"]].index.tolist()

    rows = []
    for fam in big_families:
        fam_crypt = crypt[crypt["protein_family"] == fam]["crypticity_index"]
        fam_div = diversity[diversity["protein_family"] == fam].dropna(
            subset=["crypticity_index", "n_distinct_scaffolds"])

        rho, p = (float("nan"), float("nan"))
        if len(fam_div) >= 20 and fam_div["n_distinct_scaffolds"].nunique() > 1:
            rho, p = stats.spearmanr(fam_div["crypticity_index"], fam_div["n_distinct_scaffolds"])

        rows.append({
            "protein_family": fam,
            "n_sites": len(fam_crypt),
            "median_crypticity": fam_crypt.median(),
            "mean_crypticity": fam_crypt.mean(),
            "diversity_rho": rho,
            "diversity_p": p,
        })

    fam_df = pd.DataFrame(rows).sort_values("n_sites", ascending=False)
    fam_out = f"{cfg['paths']['results']}/family_stratification.csv"
    fam_df.to_csv(fam_out, index=False)

    lines = [f"Families with >= {strat['min_family_n']} scored sites: {len(fam_df)}"]

    groups = [crypt[crypt["protein_family"] == f]["crypticity_index"].values for f in big_families]
    if len(groups) >= 2:
        h, p = stats.kruskal(*groups)
        lines.append(f"\nKruskal-Wallis, crypticity_index across all {len(groups)} families: "
                      f"H={h:.1f}, p={p:.3e} "
                      f"({'crypticity varies significantly by family' if p < 0.05 else 'no significant family effect'})")

    top_by_n = fam_df.head(strat["top_n_families"])
    lines.append(f"\nTop {strat['top_n_families']} families by n_sites (most statistically robust):")
    lines.append(top_by_n[["protein_family", "n_sites", "median_crypticity"]].to_string(index=False))

    testable = fam_df.dropna(subset=["diversity_rho"])
    top_by_rho = testable.sort_values("diversity_rho", ascending=False).head(strat["top_n_families"])
    lines.append(f"\nTop {strat['top_n_families']} families by diversity_rho "
                  f"(strongest within-family crypticity-vs-diversity signal, "
                  f"of {len(testable)} testable families):")
    lines.append(top_by_rho[["protein_family", "n_sites", "diversity_rho", "diversity_p"]].to_string(index=False))

    report = "\n".join(lines)
    print(report)
    txt_out = f"{cfg['paths']['results']}/family_stratification.txt"
    Path(txt_out).write_text(report)
    print(f"\nWrote {fam_out} + {txt_out}")


if __name__ == "__main__":
    main()
