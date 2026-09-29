"""
Collapse matched_pockets.csv (per-structure pocket observations, clustered
into cross-structure pocket_site_id) into one row per pocket site with:

    V_apo      = smallest fpocket volume among that site's apo observations
    V_liganded = largest fpocket volume among that site's holo observations
    delta_V_pocket   = V_liganded - V_apo
    crypticity_index = V_liganded / (V_apo + epsilon)

Only sites observed in BOTH at least one apo and one holo structure are
scored (crypticity is meaningless without both states to compare).

Output: results/crypticity_index.csv
    uniprot, pocket_site_id, n_apo_obs, n_holo_obs,
    v_apo_min, v_holo_max, delta_v_pocket, crypticity_index,
    holo_pdb_at_max, holo_ligand_codes_at_max
"""
import sys
import pandas as pd
from common import load_config


def main():
    cfg = load_config()
    in_path = f"{cfg['paths']['pockets']}/matched_pockets.csv"
    try:
        df = pd.read_csv(in_path)
    except FileNotFoundError:
        sys.exit(f"Missing {in_path}. Run 05_align_and_match_pockets.py first.")

    eps = cfg["crypticity"]["epsilon"]
    rows = []
    for site_id, g in df.groupby("pocket_site_id"):
        apo = g[g["is_apo"]]
        holo = g[~g["is_apo"]]
        if apo.empty or holo.empty:
            continue  # can't compute crypticity without both states at this site

        v_apo_min = apo["volume"].min()
        holo_max_row = holo.loc[holo["volume"].idxmax()]
        v_holo_max = holo_max_row["volume"]

        rows.append({
            "uniprot": g["uniprot"].iloc[0],
            "pocket_site_id": site_id,
            "n_apo_obs": len(apo),
            "n_holo_obs": len(holo),
            "v_apo_min": v_apo_min,
            "v_holo_max": v_holo_max,
            "delta_v_pocket": v_holo_max - v_apo_min,
            "crypticity_index": v_holo_max / (v_apo_min + eps),
            "druggability_apo_min_v": apo.loc[apo["volume"].idxmin(), "druggability_score"],
            "druggability_holo_max_v": holo_max_row["druggability_score"],
            "holo_pdb_at_max": holo_max_row["pdb_id"],
            "holo_ligand_codes_at_max": holo_max_row["ligand_codes"],
        })

    out = pd.DataFrame(rows).sort_values("crypticity_index", ascending=False)
    out_path = f"{cfg['paths']['results']}/crypticity_index.csv"
    out.to_csv(out_path, index=False)

    print(f"Scored {len(out)} pocket sites (with both apo and holo observations) "
          f"across {out['uniprot'].nunique()} proteins.")
    if not out.empty:
        print("Top 5 most cryptic sites:")
        print(out[["uniprot", "pocket_site_id", "v_apo_min", "v_holo_max",
                    "crypticity_index"]].head(5).to_string(index=False))
    print(f"Wrote -> {out_path}")


if __name__ == "__main__":
    main()
