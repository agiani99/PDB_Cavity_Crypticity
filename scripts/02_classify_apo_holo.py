"""
Classify each PDB entry from 01_fetch_pdb_index.py as apo (no ligand HET
codes outside the solvent/cryoprotectant blacklist) or holo (has at least
one), then keep only UniProt groups that actually have both apo and holo
members (i.e. groups where crypticity is even observable).

Output: data/index/protein_groups.csv
    uniprot, pdb_id, entity_id, resolution, is_apo, ligand_codes
"""
import sys
import pandas as pd
from common import load_config


def main():
    cfg = load_config()
    in_path = f"{cfg['paths']['index']}/pdb_groups.csv"
    try:
        df = pd.read_csv(in_path, dtype={"het_codes": str}).fillna({"het_codes": ""})
    except FileNotFoundError:
        sys.exit(f"Missing {in_path}. Run 01_fetch_pdb_index.py first.")

    blacklist = set(cfg["selection"]["solvent_blacklist"])

    def real_ligands(het_codes: str) -> str:
        codes = [c for c in het_codes.split(";") if c and c not in blacklist]
        return ";".join(sorted(set(codes)))

    df["ligand_codes"] = df["het_codes"].apply(real_ligands)
    df["is_apo"] = df["ligand_codes"] == ""

    sel = cfg["selection"]
    keep = []
    for uniprot, g in df.groupby("uniprot"):
        n_apo = g["is_apo"].sum()
        n_holo = (~g["is_apo"]).sum()
        if n_apo >= sel["min_apo_structures"] and n_holo >= sel["min_holo_structures"] \
                and len(g) >= sel["min_structures_per_protein"]:
            keep.append(uniprot)

    out = df[df["uniprot"].isin(keep)].drop(columns=["het_codes"])
    out_path = f"{cfg['paths']['index']}/protein_groups.csv"
    out.to_csv(out_path, index=False)

    print(f"{len(keep)} / {df['uniprot'].nunique()} proteins have both apo and "
          f"holo structures satisfying config.yaml [selection] thresholds.")
    print(f"Wrote {len(out)} entries -> {out_path}")


if __name__ == "__main__":
    main()
