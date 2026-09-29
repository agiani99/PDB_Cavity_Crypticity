"""
For each UniProt group: superpose every structure's CA atoms onto a
reference structure (config.yaml [pocket_matching].reference_selection),
transform every fpocket pocket centroid into that reference frame, then
greedily cluster pockets across structures into "pocket sites" — the same
physical cavity, observed (or not observed) across multiple crystal forms.

Alignment is by matching residue number (+ optional chain), which assumes
consistent PDB numbering across entries of the same protein — true for most
well-behaved entries, but a known simplification (see README "Limitations").

Output: data/pockets/matched_pockets.csv
    ...all fpocket_features.csv columns..., aligned (bool),
    aligned_centroid_x/y/z, pocket_site_id
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Superimposer
from common import load_config


def ca_coords_by_resnum(pdb_path: str) -> dict:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(Path(pdb_path).stem, pdb_path)
    coords = {}
    for chain in structure[0]:
        for residue in chain:
            if "CA" in residue:
                coords[residue.id[1]] = residue["CA"].get_coord()
    return coords


def fit_transform(ref_coords: dict, mob_coords: dict):
    """Returns (rot, tran, n_common) mapping mobile-frame coords to the
    reference frame: transformed = mobile @ rot + tran. Identity if the
    structure IS the reference."""
    common = sorted(set(ref_coords) & set(mob_coords))
    if len(common) < 3:
        return None, None, len(common)
    fixed = np.array([ref_coords[r] for r in common])
    moving = np.array([mob_coords[r] for r in common])
    sup = Superimposer()
    sup.set_atoms(_as_atoms(fixed), _as_atoms(moving))
    rot, tran = sup.rotran
    return rot, tran, len(common)


class _FakeAtom:
    """Minimal stand-in so Superimposer.set_atoms can consume raw coord arrays."""
    def __init__(self, coord):
        self._coord = np.asarray(coord, dtype=float)

    def get_coord(self):
        return self._coord

    def get_array(self):
        return self._coord

    coord = property(lambda self: self._coord)


def _as_atoms(coords):
    return [_FakeAtom(c) for c in coords]


def lining_no_chain(lining_residues: str) -> set:
    if not isinstance(lining_residues, str) or not lining_residues:
        return set()
    return {tok.split(":", 1)[1] if ":" in tok else tok for tok in lining_residues.split(";")}


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def cluster_pockets(pockets: list[dict], centroid_max_A: float, jaccard_min: float) -> list[int]:
    """Greedy single-pass clustering. Returns a cluster-id list parallel to `pockets`."""
    clusters = []  # list of {"centroid": np.array, "lining": set}
    assignment = []
    for p in pockets:
        centroid = np.array([p["aligned_centroid_x"], p["aligned_centroid_y"], p["aligned_centroid_z"]])
        lining = p["_lining_set"]
        best_idx, best_dist = None, None
        for i, c in enumerate(clusters):
            dist = np.linalg.norm(centroid - c["centroid"])
            if dist <= centroid_max_A and jaccard(lining, c["lining"]) >= jaccard_min:
                if best_dist is None or dist < best_dist:
                    best_idx, best_dist = i, dist
        if best_idx is None:
            clusters.append({"centroid": centroid, "lining": lining})
            assignment.append(len(clusters) - 1)
        else:
            assignment.append(best_idx)
    return assignment


def main():
    cfg = load_config()
    # Prefer the redundancy-reduced input (04b_deduplicate_structures.py) when
    # it exists — collapses repeat depositions of the same protein+ligand
    # complex to one representative, avoiding pseudo-replicated pocket sites.
    # Falls back to the raw fpocket output if 04b hasn't been run yet.
    deduped_pockets = Path(cfg["paths"]["pockets"]) / "fpocket_features_deduped.csv"
    deduped_idx = Path(cfg["paths"]["structures"]) / "index_deduped.csv"
    if deduped_pockets.exists() and deduped_idx.exists():
        pockets_path, idx_path = str(deduped_pockets), str(deduped_idx)
        print(f"[using deduplicated input from 04b_deduplicate_structures.py]")
    else:
        pockets_path = f"{cfg['paths']['pockets']}/fpocket_features.csv"
        idx_path = f"{cfg['paths']['structures']}/index.csv"
    try:
        pockets = pd.read_csv(pockets_path)
        idx = pd.read_csv(idx_path)
    except FileNotFoundError as e:
        sys.exit(f"Missing input ({e}). Run 04_run_fpocket.py first.")

    pm = cfg["pocket_matching"]
    all_rows = []

    for uniprot, group in idx.groupby("uniprot"):
        apo_rows = group[group["is_apo"]]
        ref_row = None
        if pm["reference_selection"] == "best_resolution_apo" and not apo_rows.empty:
            ref_row = apo_rows.sort_values("resolution").iloc[0]
        else:
            ref_row = group.sort_values("resolution").iloc[0]

        ref_coords = ca_coords_by_resnum(ref_row["clean_pdb"])
        group_pockets = pockets[pockets["uniprot"] == uniprot].copy()
        if group_pockets.empty:
            continue

        aligned_rows = []
        for pdb_id, sub in group_pockets.groupby("pdb_id"):
            struct_row = group[group["pdb_id"] == pdb_id]
            if struct_row.empty:
                continue
            mob_coords = ca_coords_by_resnum(struct_row.iloc[0]["clean_pdb"])
            rot, tran, n_common = fit_transform(ref_coords, mob_coords)
            aligned_ok = rot is not None
            for _, p in sub.iterrows():
                p = p.to_dict()
                centroid = np.array([p["centroid_x"], p["centroid_y"], p["centroid_z"]])
                if aligned_ok:
                    t_centroid = centroid @ rot + tran
                else:
                    t_centroid = centroid  # unaligned fallback: leave in original frame
                p["aligned"] = aligned_ok
                p["n_common_ca"] = n_common
                p["aligned_centroid_x"], p["aligned_centroid_y"], p["aligned_centroid_z"] = t_centroid
                p["_lining_set"] = lining_no_chain(p.get("lining_residues", ""))
                aligned_rows.append(p)

        assignment = cluster_pockets(
            aligned_rows, pm["centroid_distance_max_A"], pm["residue_jaccard_min"])
        for p, cid in zip(aligned_rows, assignment):
            p["pocket_site_id"] = f"{uniprot}_site{cid + 1}"
            del p["_lining_set"]
            all_rows.append(p)

        n_sites = len(set(assignment))
        print(f"[ok] {uniprot}: {len(group)} structures, {len(aligned_rows)} pockets "
              f"-> {n_sites} matched sites (ref={ref_row['pdb_id']})")

    out = pd.DataFrame(all_rows)
    out_path = f"{cfg['paths']['pockets']}/matched_pockets.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} pocket observations -> {out_path}")


if __name__ == "__main__":
    main()
