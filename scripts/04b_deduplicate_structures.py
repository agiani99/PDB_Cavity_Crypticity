"""
Redundancy reduction — addresses Limitation 5 from the manuscript: no
near-duplicate structure filtering was previously applied, so a protein
crystallized hundreds of times with the same bound ligand (e.g. Q6PJP8/MLI:
277 structures; P69905/HEM: 103 structures) contributed one pseudo-replicated
"observation" per re-deposition instead of one real complex.

Redundancy is defined per protein, on HOLO structures only, by their exact
sorted ligand HET-code set: two holo structures of the same protein with the
identical set of bound (non-solvent) ligands are treated as repeat
depositions of the same complex, and only the best-resolution one is kept.
Apo structures are deliberately NOT collapsed — their volume variability
across independent depositions is the actual signal hypothesis A tests, so
thinning them would remove real data, not redundancy.

This filters already-computed fpocket output (no re-downloading or
re-running fpocket) so it's cheap to re-run if the redundancy definition
needs adjusting later.

Output: data/structures/index_deduped.csv
        data/pockets/fpocket_features_deduped.csv
        results/redundancy_reduction_log.csv
"""
import sys
from pathlib import Path

import pandas as pd
from common import load_config


def main():
    cfg = load_config()
    idx_path = f"{cfg['paths']['structures']}/index.csv"
    features_path = f"{cfg['paths']['pockets']}/fpocket_features.csv"
    try:
        idx = pd.read_csv(idx_path)
        features = pd.read_csv(features_path)
    except FileNotFoundError as e:
        sys.exit(f"Missing input ({e}). Run 03_download_and_clean.py and 04_run_fpocket.py first.")

    apo = idx[idx["is_apo"]].copy()
    holo = idx[~idx["is_apo"]].copy()

    log_rows = []
    keep_pdb_ids = set(apo["pdb_id"])
    for uniprot, g in holo.groupby("uniprot"):
        for ligand_codes, dup_group in g.groupby("ligand_codes"):
            best = dup_group.sort_values("resolution").iloc[0]
            keep_pdb_ids.add(best["pdb_id"])
            n_before, n_after_this_group = len(dup_group), 1
            if n_before > 1:
                log_rows.append({
                    "uniprot": uniprot,
                    "ligand_codes": ligand_codes,
                    "n_structures_before": n_before,
                    "kept_pdb_id": best["pdb_id"],
                    "kept_resolution": best["resolution"],
                    "dropped_pdb_ids": ";".join(sorted(set(dup_group["pdb_id"]) - {best["pdb_id"]})),
                })

    idx_deduped = idx[idx["pdb_id"].isin(keep_pdb_ids)].copy()
    features_deduped = features[features["pdb_id"].isin(keep_pdb_ids)].copy()

    idx_out = f"{cfg['paths']['structures']}/index_deduped.csv"
    features_out = f"{cfg['paths']['pockets']}/fpocket_features_deduped.csv"
    log_out = f"{cfg['paths']['results']}/redundancy_reduction_log.csv"
    idx_deduped.to_csv(idx_out, index=False)
    features_deduped.to_csv(features_out, index=False)
    pd.DataFrame(log_rows).sort_values("n_structures_before", ascending=False).to_csv(log_out, index=False)

    n_holo_before, n_holo_after = len(holo), len(idx_deduped) - len(apo)
    print(f"Apo structures: {len(apo)} (kept unchanged, not deduplicated)")
    print(f"Holo structures: {n_holo_before} -> {n_holo_after} "
          f"({n_holo_before - n_holo_after} dropped as redundant, "
          f"{100 * (n_holo_before - n_holo_after) / n_holo_before:.1f}%)")
    print(f"Total structures: {len(idx)} -> {len(idx_deduped)}")
    print(f"{len(log_rows)} redundant (ligand_codes, protein) groups collapsed "
          f"(largest: see {log_out})")
    print(f"\nWrote {idx_out}\nWrote {features_out}\nWrote {log_out}")


if __name__ == "__main__":
    main()
