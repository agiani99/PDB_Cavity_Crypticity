"""
Hypothesis D (affinity), finally testable: joins each pocket site's
max-volume-holo ligand (the same PDB+HET-code pair that defines
crypticity_index in 06_compute_crypticity.py) against real experimental
ChEMBL bioactivity data via the PDBeChem bioactivity report
(config.yaml [bioactivity].report_path — an external dataset, not part of
this repo).

Deliberately uses the SAME ligand that defines crypticity_index (the one at
V_holo_max), not "the best affinity seen anywhere at this site" — keeps the
test internally consistent: one structure, one volume, one affinity.

Where a PDB+HET-code pair has multiple ChEMBL assay rows (different assays/
publications), takes the median PCHEMBL_VALUE across them (robust to a
single noisy assay) and records how many assays contributed.

Output: results/crypticity_index_with_affinity.csv
    ...crypticity_index.csv columns..., first_ligand, pchembl_value,
    n_assays, standard_type_mode
    results/affinity_correlation.txt
"""
import sys
from pathlib import Path

import pandas as pd
from scipy import stats
from common import load_config


def first_code(codes):
    codes = str(codes).split(";") if pd.notna(codes) else []
    return codes[0] if codes else None


def main():
    cfg = load_config()
    report_path = Path(cfg["bioactivity"]["report_path"])
    if not report_path.exists():
        sys.exit(f"Missing bioactivity report at {report_path} "
                  f"(check config.yaml [bioactivity].report_path).")

    try:
        crypt = pd.read_csv(f"{cfg['paths']['results']}/crypticity_index.csv")
    except FileNotFoundError:
        sys.exit("Missing results/crypticity_index.csv. Run 06_compute_crypticity.py first.")

    print(f"Loading bioactivity report from {report_path}...")
    bio = pd.read_csv(report_path, sep="\t", low_memory=False,
                       usecols=["PDB", "LIG", "PCHEMBL_VALUE", "STANDARD_TYPE"])
    bio["PDB"] = bio["PDB"].str.upper()
    bio = bio.dropna(subset=["PCHEMBL_VALUE"])

    agg = bio.groupby(["PDB", "LIG"]).agg(
        pchembl_value=("PCHEMBL_VALUE", "median"),
        n_assays=("PCHEMBL_VALUE", "size"),
        standard_type_mode=("STANDARD_TYPE", lambda s: s.mode().iat[0] if not s.mode().empty else None),
    ).reset_index()

    crypt["first_ligand"] = crypt["holo_ligand_codes_at_max"].apply(first_code)
    out = crypt.merge(agg, left_on=["holo_pdb_at_max", "first_ligand"],
                       right_on=["PDB", "LIG"], how="left").drop(columns=["PDB", "LIG"])

    out_path = f"{cfg['paths']['results']}/crypticity_index_with_affinity.csv"
    out.to_csv(out_path, index=False)

    sub = out.dropna(subset=["crypticity_index", "pchembl_value"])
    lines = [
        f"{len(sub)} / {len(out)} pocket sites have a real ChEMBL affinity value "
        f"(PCHEMBL_VALUE, median across {sub['n_assays'].sum()} total assay rows) "
        f"for their max-volume-holo ligand.",
    ]
    if len(sub) >= 10:
        rho, p = stats.spearmanr(sub["crypticity_index"], sub["pchembl_value"])
        # Statistical significance alone is a low bar at n~thousands (see hypothesis
        # B in 08_correlate_crypticity.py) — flag a significant-but-tiny effect
        # rather than call it "supported" outright.
        if p >= 0.05:
            verdict = "does not support (not significant)"
        elif abs(rho) < 0.1:
            verdict = "technically significant but negligible effect size — do not oversell this"
        else:
            verdict = "SUPPORTS"
        lines.append(f"\nD) crypticity_index vs. ligand potency (PCHEMBL_VALUE, higher = more potent):")
        lines.append(f"   n={len(sub)}, Spearman rho={rho:.3f}, p={p:.4f} ({verdict})")
    else:
        lines.append("Not enough sites with both crypticity_index and affinity to test (need >=10).")

    report = "\n".join(lines)
    print(report)
    txt_path = f"{cfg['paths']['results']}/affinity_correlation.txt"
    Path(txt_path).write_text(report)
    print(f"\nWrote {out_path} + {txt_path}")


if __name__ == "__main__":
    main()
