"""
Positive-control validation against CryptoBench (Skrhak et al. 2024,
Bioinformatics, doi:10.1093/bioinformatics/btae745) — a literature-curated
benchmark of 5,493 confirmed apo/holo cryptic-binding-site pairs (1,110
distinct UniProt accessions), the top-priority follow-up flagged in this
manuscript's Limitations §7.

(Note: the originally-suggested github.com/salilab/cryptosite repo does not
itself host a benchmark PDB list — its "data/ligands" directory is a library
of small organic solvent probe molecules used internally by CryptoSite's own
algorithm, not annotated cryptic-site structures. CryptoBench is the modern,
actively maintained, much larger successor benchmark for this exact purpose,
built specifically to evaluate methods like this one, so it was used instead.)

Two levels of validation, both against crypticity_index:

  Level 1 (protein-level, broad): is crypticity_index higher, on average,
      for proteins CryptoBench independently confirms have a cryptic site,
      vs. the rest of our dataset? Uses every protein in our dataset that
      matches a CryptoBench UniProt accession, regardless of which exact
      PDB entries CryptoBench used.

  Level 2 (site-level, precise): for the subset of CryptoBench apo PDB
      entries that are ALSO one of our own apo-structure observations,
      match CryptoBench's annotated cryptic-pocket residues (chain+resnum)
      against our matched pocket sites' lining residues at that exact
      structure (Jaccard similarity) to identify which specific
      pocket_site_id corresponds to the literature-confirmed cryptic site.
      Compares crypticity_index for these specifically-confirmed sites
      against all other scored sites.

Output: results/cryptobench_validation_summary.txt
        results/cryptobench_site_matches.csv
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd
from scipy import stats
from common import load_config

CRYPTOBENCH_PATH = "data/external/cryptobench_dataset.json"
MIN_JACCARD = 0.2  # minimum residue-set overlap to accept a site match


def parse_our_residues(lining_str: str) -> set:
    """'D:ALA196;E:ARG80' -> {('D','196'), ('E','80')}"""
    if not isinstance(lining_str, str) or not lining_str:
        return set()
    out = set()
    for tok in lining_str.split(";"):
        if ":" not in tok:
            continue
        chain, resname_resnum = tok.split(":", 1)
        m = re.search(r"(\d+)$", resname_resnum)
        if m:
            out.add((chain, m.group(1)))
    return out


def parse_cryptobench_residues(selection: list) -> set:
    """['B_12', 'B_14', ...] -> {('B','12'), ('B','14'), ...}"""
    out = set()
    for tok in selection:
        if "_" not in tok:
            continue
        chain, resnum = tok.split("_", 1)
        out.add((chain, resnum))
    return out


def jaccard(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def main():
    cfg = load_config()
    cb_path = Path(cfg["paths"]["results"]).parent / CRYPTOBENCH_PATH
    if not cb_path.exists():
        sys.exit(f"Missing {cb_path}. Download from https://osf.io/pz4a9/ "
                  f"(cryptobench/cryptobench-dataset/dataset.json).")

    try:
        crypt = pd.read_csv(f"{cfg['paths']['results']}/crypticity_index_with_family.csv")
        matched = pd.read_csv(f"{cfg['paths']['pockets']}/matched_pockets.csv")
    except FileNotFoundError as e:
        sys.exit(f"Missing input ({e}). Run 06/10/11 first.")

    cryptobench = json.loads(cb_path.read_text())
    cb_uniprots = set()
    for apo_pdb, records in cryptobench.items():
        for r in records:
            cb_uniprots.add(r["uniprot_id"])

    lines = ["=" * 70, "CryptoBench positive-control validation", "=" * 70, ""]
    lines.append(f"CryptoBench: {len(cryptobench)} apo PDB entries, "
                 f"{sum(len(v) for v in cryptobench.values())} cryptic-site records, "
                 f"{len(cb_uniprots)} distinct UniProt accessions.")

    # --- Level 1: protein-level ---
    our_uniprots = set(crypt["uniprot"].unique())
    overlap = our_uniprots & cb_uniprots
    lines.append(f"\n--- Level 1: protein-level ---")
    lines.append(f"{len(overlap)} / {len(our_uniprots)} proteins in our dataset are "
                 f"independently confirmed by CryptoBench to have a cryptic site.")

    per_protein_max_ci = crypt.groupby("uniprot")["crypticity_index"].max()
    confirmed = per_protein_max_ci[per_protein_max_ci.index.isin(overlap)]
    unconfirmed = per_protein_max_ci[~per_protein_max_ci.index.isin(overlap)]
    u, p = stats.mannwhitneyu(confirmed, unconfirmed, alternative="greater")
    lines.append(f"Per-protein max crypticity_index: CryptoBench-confirmed (n={len(confirmed)}) "
                 f"median={confirmed.median():.3f}, mean={confirmed.mean():.3f} | "
                 f"not independently confirmed (n={len(unconfirmed)}) "
                 f"median={unconfirmed.median():.3f}, mean={unconfirmed.mean():.3f}")
    lines.append(f"Mann-Whitney U (H1: confirmed > unconfirmed): p={p:.3e} "
                 f"({'SUPPORTS' if p < 0.05 else 'does not support'} — crypticity_index runs "
                 f"higher for proteins with an independently literature-confirmed cryptic site)")

    # --- Level 2: site-level, exact residue matching ---
    our_apo = matched[matched["is_apo"]].copy()
    our_apo["pdb_upper"] = our_apo["pdb_id"].str.upper()
    lines.append(f"\n--- Level 2: site-level (exact apo-structure residue matching) ---")

    site_rows = []
    for apo_pdb, records in cryptobench.items():
        candidates = our_apo[our_apo["pdb_upper"] == apo_pdb.upper()]
        if candidates.empty:
            continue
        for r in records:
            cb_res = parse_cryptobench_residues(r["apo_pocket_selection"])
            best_site, best_j = None, 0.0
            for _, row in candidates.iterrows():
                our_res = parse_our_residues(row["lining_residues"])
                j = jaccard(cb_res, our_res)
                if j > best_j:
                    best_j, best_site = j, row["pocket_site_id"]
            if best_site is not None and best_j >= MIN_JACCARD:
                site_rows.append({
                    "apo_pdb": apo_pdb, "uniprot": r["uniprot_id"],
                    "holo_pdb": r["holo_pdb_id"], "pocket_site_id": best_site,
                    "jaccard": best_j,
                })

    site_df = pd.DataFrame(site_rows).drop_duplicates(subset=["pocket_site_id"])
    lines.append(f"{len(cryptobench)} CryptoBench apo entries checked against our apo structures; "
                 f"{our_apo['pdb_upper'].isin([p.upper() for p in cryptobench]).sum()} of our apo "
                 f"observations are an exact CryptoBench apo entry; "
                 f"{len(site_df)} distinct pocket_site_id's matched a CryptoBench-annotated cryptic "
                 f"pocket at Jaccard >= {MIN_JACCARD}.")

    if len(site_df) >= 5:
        site_df = site_df.merge(crypt[["pocket_site_id", "crypticity_index"]],
                                 on="pocket_site_id", how="left").dropna(subset=["crypticity_index"])
        confirmed_sites = site_df["crypticity_index"]
        other_sites = crypt.loc[~crypt["pocket_site_id"].isin(site_df["pocket_site_id"]), "crypticity_index"]
        u2, p2 = stats.mannwhitneyu(confirmed_sites, other_sites, alternative="greater")
        lines.append(f"\nCrypticity index at CryptoBench-confirmed sites (n={len(confirmed_sites)}): "
                     f"median={confirmed_sites.median():.3f}, mean={confirmed_sites.mean():.3f}")
        lines.append(f"Crypticity index at all other scored sites (n={len(other_sites)}): "
                     f"median={other_sites.median():.3f}, mean={other_sites.mean():.3f}")
        lines.append(f"Mann-Whitney U (H1: confirmed sites score higher): p={p2:.3e} "
                     f"({'SUPPORTS' if p2 < 0.05 else 'does not support'})")
    else:
        lines.append("Too few site-level matches for a meaningful test.")

    report = "\n".join(lines)
    print(report)
    out_txt = f"{cfg['paths']['results']}/cryptobench_validation_summary.txt"
    Path(out_txt).write_text(report, encoding="utf-8")
    site_df.to_csv(f"{cfg['paths']['results']}/cryptobench_site_matches.csv", index=False)
    print(f"\nWrote {out_txt} + results/cryptobench_site_matches.csv")


if __name__ == "__main__":
    main()
