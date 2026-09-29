"""
Download each PDB entry in data/index/protein_groups.csv and strip it down
to protein-only coordinates (waters, ions, and the bound ligand itself all
removed) so fpocket sees a comparable apo-style cavity for every structure,
whether it was solved apo or holo.

Output: data/structures/<pdb_id>_clean.pdb  (one per entry)
        data/structures/index.csv  (pdb_id, uniprot, is_apo, resolution, clean_pdb)
"""
import sys
import pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from common import load_config, download_pdb_file


class ProteinOnly(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "  # standard ATOM records only, no HETATM/water


def clean_structure(raw_path: Path, out_pdb: Path):
    parser = MMCIFParser(QUIET=True) if raw_path.suffix == ".cif" else PDBParser(QUIET=True)
    structure = parser.get_structure(raw_path.stem, str(raw_path))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_pdb), ProteinOnly())


def main():
    cfg = load_config()
    in_path = f"{cfg['paths']['index']}/protein_groups.csv"
    try:
        df = pd.read_csv(in_path)
    except FileNotFoundError:
        sys.exit(f"Missing {in_path}. Run 02_classify_apo_holo.py first.")

    raw_dir = Path(cfg["paths"]["structures"]) / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    clean_dir = Path(cfg["paths"]["structures"])

    rows = []
    for _, row in df.drop_duplicates(subset=["pdb_id"]).iterrows():
        pdb_id = row["pdb_id"]
        try:
            raw_path = download_pdb_file(pdb_id, str(raw_dir), cfg)
        except Exception as e:
            print(f"[skip] {pdb_id}: download failed ({e})")
            continue

        clean_path = clean_dir / f"{pdb_id}_clean.pdb"
        try:
            clean_structure(raw_path, clean_path)
        except Exception as e:
            print(f"[skip] {pdb_id}: cleaning failed ({e})")
            continue

        rows.append({
            "pdb_id": pdb_id,
            "uniprot": row["uniprot"],
            "is_apo": row["is_apo"],
            "ligand_codes": row["ligand_codes"],
            "resolution": row["resolution"],
            "clean_pdb": str(clean_path),
        })
        print(f"[ok] {pdb_id} ({'apo' if row['is_apo'] else 'holo'}) -> {clean_path.name}")

    out = pd.DataFrame(rows)
    out_path = clean_dir / "index.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} cleaned structures -> {out_path}")


if __name__ == "__main__":
    main()
