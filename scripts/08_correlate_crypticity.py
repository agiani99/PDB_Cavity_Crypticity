"""
Test whether structural variability predicts cryptic-pocket formation and
its downstream ligand properties, across the four axes from the original
question:

  A) ligand presence   — do sites whose APO volume varies more across
                          crystal forms end up binding a ligand more often?
                          (Mann-Whitney U + logistic regression)
  B) ligand size        — does crypticity_index correlate with the bound
                          ligand's heavy-atom count / MW? (Spearman)
  C) chemical diversity — do more-cryptic sites accommodate a wider variety
                          of scaffolds across their holo structures? (Spearman)
  D) affinity            — NOT COMPUTED: no experimental Kd/IC50/Ki source is
                          wired up yet (most PDB ligands have none anyway).
                          Wire in PDBbind/BindingMOAD and re-run this stage
                          if you need this arm.

Output: results/correlation_summary.txt, results/site_level_analysis.csv
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
from common import load_config


def hypothesis_ligand_presence(matched: pd.DataFrame) -> tuple[str, pd.DataFrame]:
    lines = ["--- A) Ligand presence vs. apo structural variability ---"]
    per_site = []
    for site_id, g in matched.groupby("pocket_site_id"):
        apo = g[g["is_apo"]]
        has_ligand = (~g["is_apo"]).any()
        if len(apo) < 2:
            continue  # need >=2 apo observations to see any variability at all
        per_site.append({
            "pocket_site_id": site_id,
            "uniprot": g["uniprot"].iloc[0],
            "n_apo_obs": len(apo),
            "apo_volume_range": apo["volume"].max() - apo["volume"].min(),
            "apo_volume_std": apo["volume"].std(),
            "has_ligand": bool(has_ligand),
        })
    df = pd.DataFrame(per_site)
    if df.empty or df["has_ligand"].nunique() < 2:
        lines.append("Not enough sites with >=2 apo observations in both outcome "
                      "classes to test this arm yet.")
        return "\n".join(lines), df

    with_l = df[df["has_ligand"]]["apo_volume_range"].dropna()
    without_l = df[~df["has_ligand"]]["apo_volume_range"].dropna()
    if len(with_l) >= 3 and len(without_l) >= 3:
        u, p = stats.mannwhitneyu(with_l, without_l, alternative="greater")
        lines.append(f"n(ligand-bound sites)={len(with_l)}, n(apo-only sites)={len(without_l)}")
        lines.append(f"Mann-Whitney U (H1: apo volume range is larger for sites that "
                      f"go on to bind a ligand): U={u:.1f}, p={p:.4f} "
                      f"({'SUPPORTS' if p < 0.05 else 'does not support'} the hypothesis)")
        lines.append(f"median apo_volume_range: ligand-bound={with_l.median():.1f} A^3, "
                      f"apo-only={without_l.median():.1f} A^3")
    else:
        lines.append("Not enough sites per group (need >=3 each) for Mann-Whitney U.")

    try:
        import statsmodels.api as sm
        X = sm.add_constant(df["apo_volume_range"].fillna(0.0))
        model = sm.Logit(df["has_ligand"].astype(int), X).fit(disp=0)
        lines.append("\nLogistic regression: has_ligand ~ apo_volume_range")
        lines.append(str(model.summary2().tables[1]))
    except Exception as e:
        lines.append(f"[skip] logistic regression failed: {e}")

    return "\n".join(lines), df


def hypothesis_size_and_diversity(crypticity: pd.DataFrame, ligands: pd.DataFrame,
                                    matched: pd.DataFrame) -> str:
    lines = ["\n--- B) Ligand size vs. crypticity_index ---"]
    lig_lookup = ligands.set_index("het_code")

    def first_code(codes):
        codes = str(codes).split(";") if pd.notna(codes) else []
        return codes[0] if codes else None

    crypt = crypticity.copy()
    crypt["first_ligand"] = crypt["holo_ligand_codes_at_max"].apply(first_code)
    crypt = crypt.merge(lig_lookup[["heavy_atoms", "mol_wt", "clogp", "tpsa"]],
                         left_on="first_ligand", right_index=True, how="left")

    for col, label in [("heavy_atoms", "heavy-atom count"), ("mol_wt", "molecular weight")]:
        sub = crypt.dropna(subset=["crypticity_index", col])
        if len(sub) >= 5:
            rho, p = stats.spearmanr(sub["crypticity_index"], sub[col])
            lines.append(f"crypticity_index vs {label}: n={len(sub)}, "
                         f"Spearman rho={rho:.3f}, p={p:.4f} "
                         f"({'SUPPORTS' if p < 0.05 else 'does not support'})")
        else:
            lines.append(f"crypticity_index vs {label}: not enough data (n={len(sub)}).")

    lines.append("\n--- C) Chemical diversity of bound ligands vs. crypticity_index ---")
    het_to_scaffold = lig_lookup["murcko_scaffold_smiles"].to_dict()
    holo = matched[~matched["is_apo"]]
    diversity_rows = []
    for site_id, g in holo.groupby("pocket_site_id"):
        codes = set()
        for c in g["ligand_codes"].dropna():
            codes.update(str(c).split(";"))
        codes.discard("")
        scaffolds = {het_to_scaffold.get(c) for c in codes} - {None}
        diversity_rows.append({
            "pocket_site_id": site_id,
            "n_distinct_ligands": len(codes),
            "n_distinct_scaffolds": len(scaffolds),
        })
    div = pd.DataFrame(diversity_rows).merge(
        crypticity[["pocket_site_id", "crypticity_index"]], on="pocket_site_id", how="inner")
    sub = div[div["n_distinct_ligands"] >= 1]
    if len(sub) >= 5 and sub["n_distinct_scaffolds"].nunique() > 1:
        rho, p = stats.spearmanr(sub["crypticity_index"], sub["n_distinct_scaffolds"])
        lines.append(f"crypticity_index vs n_distinct_scaffolds: n={len(sub)}, "
                     f"Spearman rho={rho:.3f}, p={p:.4f} "
                     f"({'SUPPORTS' if p < 0.05 else 'does not support'})")
    else:
        lines.append(f"Not enough sites with multiple holo ligands to test diversity "
                     f"(n={len(sub)}, most sites likely have only 1 distinct ligand). "
                     f"This arm needs a bigger/broader structure pull to have power.")

    lines.append("\n--- D) Binding affinity vs. crypticity_index ---")
    lines.append("NOT COMPUTED — no experimental affinity source wired up (see docstring).")

    return "\n".join(lines), div


def main():
    cfg = load_config()
    try:
        matched = pd.read_csv(f"{cfg['paths']['pockets']}/matched_pockets.csv")
        crypticity = pd.read_csv(f"{cfg['paths']['results']}/crypticity_index.csv")
        ligands = pd.read_csv(f"{cfg['paths']['ligands']}/ligand_descriptors.csv")
    except FileNotFoundError as e:
        sys.exit(f"Missing input ({e}). Run 05/06/07 first.")

    report_a, site_df = hypothesis_ligand_presence(matched)
    report_bc, div_df = hypothesis_size_and_diversity(crypticity, ligands, matched)

    full_report = report_a + "\n" + report_bc + "\n"
    print(full_report)

    out_txt = f"{cfg['paths']['results']}/correlation_summary.txt"
    Path(out_txt).write_text(full_report)

    site_df.to_csv(f"{cfg['paths']['results']}/site_level_analysis.csv", index=False)
    div_df.to_csv(f"{cfg['paths']['results']}/diversity_analysis.csv", index=False)
    print(f"\nWrote {out_txt} + site_level_analysis.csv + diversity_analysis.csv")


if __name__ == "__main__":
    main()
