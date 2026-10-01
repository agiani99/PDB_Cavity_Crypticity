"""
Second, independent positive-control validation — this time against Lacuna
(Moore 2026, bioRxiv 10.64898/2026.08.14.744956), a cryptic-pocket detector
that takes a *single* static structure and generates its own conformational
ensemble (elastic-network normal-mode analysis by default, no GPU/force
field needed), detects pockets per conformer, and reports a continuous
per-site crypticity score plus a binary cryptic flag.

This is a methodologically orthogonal check to CryptoBench (15_...py):
CryptoBench is literature ground truth; Lacuna is a *second computational
method* that infers crypticity from a single apo structure's own predicted
flexibility, never seeing our multi-structure PDB variability data at all.
Agreement between the two independent computational approaches is evidence
that crypticity_index reflects a real structural property of these apo
structures, not an artifact specific to our own pipeline's volume-variability
calculation.

Stage 1 (run): for each protein's reference apo structure (the lowest-
resolution apo structure, matching the reference choice already used in
05_align_and_match_pockets.py), run `lacuna discover` (NMA backend, no GPU)
and cache the resulting pocket_report.json. Resumable — already-run
proteins are skipped.

Stage 2 (match + correlate): match each Lacuna-detected pocket's lining
residues against our own matched pockets at that same reference structure
(same chain+resnum Jaccard approach as 15_validate_against_cryptobench.py),
then correlate Lacuna's crypticity score against crypticity_index for
matched sites, and compare crypticity_index between sites Lacuna flags as
cryptic vs. not.

Output: results/lacuna_validation_summary.txt
        results/lacuna_site_matches.csv
"""
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from scipy import stats
from common import load_config

LACUNA_REPORTS_DIR = "data/external/lacuna_reports"
N_CONFORMERS = 15
MIN_JACCARD = 0.2


def parse_our_residues(lining_str: str) -> set:
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


def parse_lacuna_residues(lining_list: list) -> set:
    """'ALA65:A' -> {('A', '65')}"""
    out = set()
    for tok in lining_list:
        if ":" not in tok:
            continue
        resname_resnum, chain = tok.rsplit(":", 1)
        m = re.search(r"(\d+)$", resname_resnum)
        if m:
            out.add((chain, m.group(1)))
    return out


def jaccard(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def run_lacuna(cfg, reference_apo: pd.DataFrame, reports_dir: Path):
    reports_dir.mkdir(parents=True, exist_ok=True)
    n_ok, n_fail, n_skip = 0, 0, 0
    for _, row in reference_apo.iterrows():
        out_path = reports_dir / f"{row['uniprot']}.json"
        if out_path.exists():
            n_skip += 1
            continue
        try:
            result = subprocess.run(
                ["lacuna", "discover", row["clean_pdb"], "-b", "nma",
                 "-n", str(N_CONFORMERS), "--quiet", "-o", "lacuna_tmp_out"],
                capture_output=True, text=True, timeout=120, cwd=reports_dir.parent,
            )
            tmp_report = reports_dir.parent / "lacuna_tmp_out" / "pocket_report.json"
            if result.returncode != 0 or not tmp_report.exists():
                print(f"[fail] {row['uniprot']} ({row['pdb_id']}): {result.stderr[-300:]}")
                n_fail += 1
                continue
            report = json.loads(tmp_report.read_text())
            report["_uniprot"] = row["uniprot"]
            report["_ref_pdb_id"] = row["pdb_id"]
            out_path.write_text(json.dumps(report), encoding="utf-8")
            n_ok += 1
        except Exception as e:
            print(f"[fail] {row['uniprot']} ({row['pdb_id']}): {e}")
            n_fail += 1
        if (n_ok + n_fail) % 25 == 0 and (n_ok + n_fail) > 0:
            print(f"  ...{n_ok} ok, {n_fail} failed, {n_skip} skipped so far")
    print(f"Lacuna runs: {n_ok} ok, {n_fail} failed, {n_skip} already cached.")


def main():
    cfg = load_config()
    root = Path(cfg["paths"]["results"]).parent
    reports_dir = root / LACUNA_REPORTS_DIR

    try:
        crypt = pd.read_csv(f"{cfg['paths']['results']}/crypticity_index_with_family.csv")
        matched = pd.read_csv(f"{cfg['paths']['pockets']}/matched_pockets.csv")
        idx = pd.read_csv(f"{cfg['paths']['structures']}/index_deduped.csv")
    except FileNotFoundError as e:
        sys.exit(f"Missing input ({e}). Run 04b/06/10 first.")

    scored_proteins = set(crypt["uniprot"].unique())
    apo = idx[idx["is_apo"] & idx["uniprot"].isin(scored_proteins)].copy()
    reference_apo = apo.sort_values("resolution").groupby("uniprot").first().reset_index()
    print(f"{len(reference_apo)} proteins with a reference apo structure and >=1 scored site.")

    if "--run" in sys.argv or not any(reports_dir.glob("*.json")):
        run_lacuna(cfg, reference_apo, reports_dir)
    else:
        print(f"Using {len(list(reports_dir.glob('*.json')))} cached Lacuna reports "
              f"(pass --run to fill in any missing ones).")
        run_lacuna(cfg, reference_apo, reports_dir)

    # --- match Lacuna pockets to our pocket_site_id's ---
    our_apo = matched[matched["is_apo"]].copy()
    site_rows = []
    n_reports = 0
    for report_path in reports_dir.glob("*.json"):
        report = json.loads(report_path.read_text())
        n_reports += 1
        uniprot, ref_pdb = report["_uniprot"], report["_ref_pdb_id"]
        candidates = our_apo[our_apo["pdb_id"].str.upper() == ref_pdb.upper()]
        if candidates.empty:
            continue
        for pocket in report["pockets"]:
            lac_res = parse_lacuna_residues(pocket["lining_residues"])
            best_site, best_j = None, 0.0
            for _, row in candidates.iterrows():
                our_res = parse_our_residues(row["lining_residues"])
                j = jaccard(lac_res, our_res)
                if j > best_j:
                    best_j, best_site = j, row["pocket_site_id"]
            if best_site is not None and best_j >= MIN_JACCARD:
                site_rows.append({
                    "uniprot": uniprot, "ref_pdb": ref_pdb,
                    "pocket_site_id": best_site, "jaccard": best_j,
                    "lacuna_rank": pocket["rank"], "lacuna_crypticity": pocket["crypticity"],
                    "lacuna_cryptic_flag": pocket["cryptic"],
                    "lacuna_druggability": pocket["druggability"],
                })

    site_df = pd.DataFrame(site_rows)
    if not site_df.empty:
        site_df = site_df.sort_values("jaccard", ascending=False).drop_duplicates(subset=["pocket_site_id"])
        site_df = site_df.merge(crypt[["pocket_site_id", "crypticity_index"]],
                                 on="pocket_site_id", how="left").dropna(subset=["crypticity_index"])

    lines = ["=" * 70, "Lacuna positive-control validation (independent computational method)",
             "=" * 70, "",
             f"Lacuna reports available: {n_reports} proteins ({N_CONFORMERS} NMA conformers each).",
             f"Matched {len(site_df)} of our pocket_site_id's to a Lacuna-detected pocket "
             f"at the same reference apo structure (Jaccard >= {MIN_JACCARD})."]

    if len(site_df) >= 10:
        rho, p = stats.spearmanr(site_df["crypticity_index"], site_df["lacuna_crypticity"])
        lines.append(f"\ncrypticity_index vs. Lacuna's own crypticity score: "
                     f"n={len(site_df)}, Spearman rho={rho:.3f}, p={p:.3e} "
                     f"({'SUPPORTS' if p < 0.05 and rho > 0 else 'does not support'} agreement "
                     f"between the two independent methods)")

        lac_cryptic = site_df.loc[site_df["lacuna_cryptic_flag"], "crypticity_index"]
        lac_not = site_df.loc[~site_df["lacuna_cryptic_flag"], "crypticity_index"]
        if len(lac_cryptic) >= 5 and len(lac_not) >= 5:
            u, p2 = stats.mannwhitneyu(lac_cryptic, lac_not, alternative="greater")
            lines.append(f"\nOur crypticity_index where Lacuna's own binary flag says cryptic "
                         f"(n={len(lac_cryptic)}): median={lac_cryptic.median():.3f} vs. "
                         f"not flagged (n={len(lac_not)}): median={lac_not.median():.3f}")
            lines.append(f"Mann-Whitney U (H1: Lacuna-flagged sites score higher on our metric): "
                         f"p={p2:.3e} ({'SUPPORTS' if p2 < 0.05 else 'does not support'})")
    else:
        lines.append("Too few matches for a meaningful correlation test.")

    report_text = "\n".join(lines)
    print(report_text)
    Path(f"{cfg['paths']['results']}/lacuna_validation_summary.txt").write_text(report_text, encoding="utf-8")
    site_df.to_csv(f"{cfg['paths']['results']}/lacuna_site_matches.csv", index=False)
    print(f"\nWrote results/lacuna_validation_summary.txt + results/lacuna_site_matches.csv")


if __name__ == "__main__":
    main()
