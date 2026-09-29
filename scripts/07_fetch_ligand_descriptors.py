"""
For every distinct ligand HET code seen in matched_pockets.csv's holo rows,
fetch its formula/weight/SMILES from the RCSB chemical component dictionary
and compute RDKit descriptors (heavy atom count, MW, logP, TPSA, Murcko
scaffold — the scaffold is the chemical-diversity handle used downstream).

A small fraction of ligands (~0.5%, mostly organometallic sandwich complexes
and boron-cluster carboranes) fail RDKit's default sanitization outright.
These are deliberately NOT force-fitted with an exotic parser — heavy_atoms
and mol_wt are recovered where possible via a permissive (sanitize=False +
UpdatePropertyCache) fallback, since atom count/weight don't depend on
valence or aromaticity perception, but clogp/tpsa/murcko_scaffold_smiles
(which do) are left null rather than guessed at. `descriptor_quality` marks
which path each row took, so the ~0.5% partial rows stay identifiable
downstream instead of silently blending in.

Output: data/ligands/ligand_descriptors.csv
    het_code, formula, formula_weight, smiles, heavy_atoms, mol_wt,
    clogp, tpsa, murcko_scaffold_smiles, descriptor_quality
"""
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from rdkit.Chem.Scaffolds import MurckoScaffold
from common import load_config, rcsb_chem_comp_descriptors


def parse_mol_permissive(smiles):
    """Returns (mol, quality). Tries a fully-sanitized parse first; falls
    back to sanitize=False (still gives atom count/weight, not much else)
    rather than dropping the ligand outright."""
    if not smiles:
        return None, "no_smiles"
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol, "full"
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None, "unparseable"
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    return mol, "partial"


def main():
    cfg = load_config()
    in_path = f"{cfg['paths']['pockets']}/matched_pockets.csv"
    try:
        df = pd.read_csv(in_path)
    except FileNotFoundError:
        sys.exit(f"Missing {in_path}. Run 05_align_and_match_pockets.py first.")

    holo = df[~df["is_apo"]]
    het_codes = set()
    for codes in holo["ligand_codes"].dropna():
        het_codes.update(c for c in str(codes).split(";") if c)

    if not het_codes:
        sys.exit("No ligand HET codes found among holo structures.")

    print(f"Fetching chemical component data for {len(het_codes)} distinct ligands...")
    comp_info = rcsb_chem_comp_descriptors(sorted(het_codes), cfg)

    rows = []
    quality_counts = {}
    for code, info in comp_info.items():
        row = {"het_code": code, **info}
        mol, quality = parse_mol_permissive(info.get("smiles"))
        quality_counts[quality] = quality_counts.get(quality, 0) + 1
        row["descriptor_quality"] = quality
        row.update({"heavy_atoms": None, "mol_wt": None, "clogp": None,
                    "tpsa": None, "murcko_scaffold_smiles": None})
        if mol is not None:
            row["heavy_atoms"] = mol.GetNumHeavyAtoms()
            try:
                row["mol_wt"] = Descriptors.MolWt(mol)
            except Exception:
                pass
        if quality == "full":
            row["clogp"] = Crippen.MolLogP(mol)
            row["tpsa"] = Descriptors.TPSA(mol)
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                row["murcko_scaffold_smiles"] = Chem.MolToSmiles(scaffold)
            except Exception:
                pass
        rows.append(row)

    out = pd.DataFrame(rows)
    out_path = f"{cfg['paths']['ligands']}/ligand_descriptors.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} ligand descriptor rows -> {out_path}")
    print(f"Descriptor quality breakdown: {quality_counts}")


if __name__ == "__main__":
    main()
