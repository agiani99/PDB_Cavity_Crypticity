"""
Hypothesis E: apo-state druggability score predicts crypticity — the most
direct test of this project's own motivating premise (a standard pocket
detector calls a cryptic site "weak/non-druggable" in its unliganded form,
yet a real ligand later proves it wrong). fpocket's own 0-1 druggability
score was already computed per matched site in 06_compute_crypticity.py
(druggability_apo_min_v / druggability_holo_max_v) but never analyzed.

Three angles, all against crypticity_index:
  E1) Spearman correlation: crypticity_index vs. apo-state druggability
      score (pooled). Predicts negative: more cryptic -> lower apo score.
  E2) Mann-Whitney U: apo druggability score, highly cryptic (CI>=5) vs.
      minimally cryptic (CI<2) sites — using this manuscript's own §3.4
      narrative thresholds, so the comparison lines up with how "cryptic"
      is already defined elsewhere in the paper.
  E3) Spearman correlation: crypticity_index vs. delta-druggability
      (holo_max - apo_min). This is closer to an internal-consistency
      check than a new independent finding, since fpocket's own
      druggability score is itself partly a function of pocket volume —
      reported for completeness, not oversold as novel.

Caveat carried through to the writeup: druggability_apo_min_v /
druggability_holo_max_v are fpocket's score for the SPECIFIC matched pocket
at that site, not necessarily the top-ranked pocket in that structure —
most individual sites (most are minor surface pockets, not the one "real"
binding site fpocket would rank #1) score near 0 on this scale regardless
of crypticity, so this is a real but small effect, not a dramatic one.

Output: results/druggability_summary.txt, results/druggability_analysis.csv
"""
import sys
from pathlib import Path

import pandas as pd
from scipy import stats
from common import load_config


def fmt_p(p: float) -> str:
    """scipy underflows p to exactly 0.0 for very strong/large-n correlations."""
    return f"p={p:.3e}" if p > 0 else "p<1e-300"


def main():
    cfg = load_config()
    try:
        c = pd.read_csv(f"{cfg['paths']['results']}/crypticity_index.csv")
    except FileNotFoundError:
        sys.exit("Missing results/crypticity_index.csv. Run 06_compute_crypticity.py first.")

    lines = []

    # E1: pooled correlation, apo druggability vs crypticity
    sub = c.dropna(subset=["druggability_apo_min_v", "crypticity_index"])
    rho, p = stats.spearmanr(sub["druggability_apo_min_v"], sub["crypticity_index"])
    lines.append("--- E1) Apo-state druggability score vs. crypticity_index (pooled) ---")
    lines.append(f"n={len(sub)}, Spearman rho={rho:.3f}, p={p:.3e} "
                 f"({'SUPPORTS' if p < 0.05 and rho < 0 else 'does not support'} "
                 f"the 'detector underrates cryptic sites' hypothesis)")
    lines.append(f"(median apo druggability score across all sites: {sub['druggability_apo_min_v'].median():.4f} "
                 f"— heavily zero-inflated: most matched sites are minor surface pockets, "
                 f"not the single top-ranked pocket in their structure)")

    # E2: highly cryptic (CI>=5) vs minimally cryptic (CI<2) -- this manuscript's own thresholds
    cryptic = c[c["crypticity_index"] >= 5]["druggability_apo_min_v"].dropna()
    noncryptic = c[c["crypticity_index"] < 2]["druggability_apo_min_v"].dropna()
    u, p2 = stats.mannwhitneyu(noncryptic, cryptic, alternative="greater")
    lines.append("\n--- E2) Apo druggability: highly cryptic (CI>=5) vs. minimally cryptic (CI<2) ---")
    lines.append(f"n(highly cryptic)={len(cryptic)}, n(minimally cryptic)={len(noncryptic)}")
    lines.append(f"median: highly cryptic={cryptic.median():.4f}, minimally cryptic={noncryptic.median():.4f}")
    lines.append(f"mean:   highly cryptic={cryptic.mean():.4f}, minimally cryptic={noncryptic.mean():.4f}")
    lines.append(f"Mann-Whitney U (H1: minimally-cryptic sites score higher): U={u:.1f}, p={p2:.3e} "
                 f"({'SUPPORTS' if p2 < 0.05 else 'does not support'})")

    # E3: crypticity vs change in druggability score (internal consistency, not independent)
    c["delta_druggability"] = c["druggability_holo_max_v"] - c["druggability_apo_min_v"]
    sub3 = c.dropna(subset=["delta_druggability", "crypticity_index"])
    rho3, p3 = stats.spearmanr(sub3["crypticity_index"], sub3["delta_druggability"])
    lines.append("\n--- E3) crypticity_index vs. delta-druggability (holo_max - apo_min) ---")
    lines.append(f"n={len(sub3)}, Spearman rho={rho3:.3f}, {fmt_p(p3)} "
                 f"(NOTE: not independent evidence — fpocket's own druggability score is "
                 f"partly volume-derived, so this is closer to an internal-consistency check "
                 f"than a new finding)")

    report = "\n".join(lines)
    print(report)
    txt_path = f"{cfg['paths']['results']}/druggability_summary.txt"
    Path(txt_path).write_text(report)

    out_path = f"{cfg['paths']['results']}/druggability_analysis.csv"
    c[["uniprot", "pocket_site_id", "crypticity_index", "druggability_apo_min_v",
       "druggability_holo_max_v", "delta_druggability"]].to_csv(out_path, index=False)
    print(f"\nWrote {txt_path} + {out_path}")


if __name__ == "__main__":
    main()
