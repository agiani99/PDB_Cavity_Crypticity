"""
Run fpocket on every cleaned structure from 03_download_and_clean.py and
parse every detected pocket (volume, druggability, centroid, lining
residues) — unlike the hot-spot mimicry pipeline, there's no single known
site to pick here; pocket correspondence across structures is resolved
later in 05_align_and_match_pockets.py.

Output: data/pockets/fpocket_features.csv
    pdb_id, uniprot, is_apo, pocket_id, druggability_score, volume,
    hydrophobicity_score, polarity_score, centroid_x/y/z, lining_residues
"""
import re
import subprocess
import sys
from pathlib import Path
import pandas as pd
from common import load_config, run_fpocket, fpocket_output_dir


def parse_pocket_info(info_path: Path) -> list[dict]:
    pockets = []
    if not info_path.exists():
        return pockets
    text = info_path.read_text()
    blocks = re.split(r"Pocket\s+(\d+)\s*:\n", text)[1:]
    for pocket_id, body in zip(blocks[0::2], blocks[1::2]):
        fields = {}
        for line in body.strip().splitlines():
            if ":" not in line:
                continue
            key, val = line.split(":", 1)
            fields[key.strip()] = val.strip()
        pockets.append({
            "pocket_id": int(pocket_id),
            "druggability_score": float(fields.get("Druggability Score", "nan")),
            "volume": float(fields.get("Volume", "nan")),
            "hydrophobicity_score": float(fields.get("Hydrophobicity score", "nan")),
            "polarity_score": float(fields.get("Polarity score", "nan")),
        })
    return pockets


def pocket_geometry(pocket_atm_pdb: Path):
    """Centroid of the pocket's alpha-sphere-adjacent atoms, and the set of
    lining residues as (chain, resnum, resname) triples."""
    xs, ys, zs, residues = [], [], [], set()
    if not pocket_atm_pdb.exists():
        return None, set()
    for line in pocket_atm_pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
            chain = line[21].strip()
            resnum = line[22:26].strip()
            resname = line[17:20].strip()
            residues.add(f"{chain}:{resname}{resnum}")
    if not xs:
        return None, set()
    centroid = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
    return centroid, residues


def main():
    cfg = load_config()
    idx_path = f"{cfg['paths']['structures']}/index.csv"
    try:
        idx = pd.read_csv(idx_path)
    except FileNotFoundError:
        sys.exit(f"Missing {idx_path}. Run 03_download_and_clean.py first.")

    binary = cfg["fpocket"]["binary"]
    rows = []
    lining_lookup = {}  # (pdb_id, pocket_id) -> set of lining residue strings

    for _, row in idx.iterrows():
        pdb_path = row["clean_pdb"]
        try:
            run_fpocket(binary, pdb_path)
            out_dir = fpocket_output_dir(pdb_path)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"[skip] fpocket failed on {pdb_path}: {e}")
            continue

        stem = Path(pdb_path).stem
        info_path = out_dir / f"{stem}_info.txt"
        pockets = parse_pocket_info(info_path)
        if not pockets:
            print(f"[skip] no pockets parsed for {pdb_path}")
            continue

        pockets_dir = out_dir / "pockets"
        for p in pockets:
            atm_pdb = pockets_dir / f"pocket{p['pocket_id']}_atm.pdb"
            centroid, residues = pocket_geometry(atm_pdb)
            if centroid is None:
                continue
            lining_lookup[(row["pdb_id"], p["pocket_id"])] = residues
            rows.append({
                "pdb_id": row["pdb_id"],
                "uniprot": row["uniprot"],
                "is_apo": row["is_apo"],
                "ligand_codes": row["ligand_codes"],
                "resolution": row["resolution"],
                **p,
                "centroid_x": round(centroid[0], 3),
                "centroid_y": round(centroid[1], 3),
                "centroid_z": round(centroid[2], 3),
                "lining_residues": ";".join(sorted(residues)),
            })
        print(f"[ok] {row['pdb_id']}: {len(pockets)} pockets")

    out = pd.DataFrame(rows)
    out_path = f"{cfg['paths']['pockets']}/fpocket_features.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} pockets across {out['pdb_id'].nunique()} structures -> {out_path}")


if __name__ == "__main__":
    main()
