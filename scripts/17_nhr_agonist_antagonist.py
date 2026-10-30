"""
Hypothesis F: pocket plasticity of nuclear hormone receptors (NHRs) depends on ligand efficacy.

Textbook expectation (e.g. ER-alpha, PPAR-gamma, ROR-gamma): agonists lock helix 12 over the
ligand-binding pocket (compact, stable pocket); antagonists / inverse agonists displace or
destabilise helix 12 (different, looser pocket). We ask three questions at matched pocket sites
that are populated by BOTH agonist-bound and antagonist-like (antagonist + inverse agonist)
holo structures:

  F1  Pocket volume: is it different between agonist- and antagonist-bound structures?
  F2  Plasticity: is the structure-to-structure spread of the volume (robust relative MAD,
      equal-sample-size resampling) smaller in agonist-bound structures?
  F3  Crypticity: is the ligand-class-specific crypticity index
      CI_class = median(V_holo, class) / (V_apo_min + eps) different between classes?

Efficacy labels (the weak point -- PDB has no efficacy field):
  1. text mining of the RCSB entry title + keywords (agonist / partial agonist / antagonist /
     inverse agonist; entries naming both agonist and antagonist are discarded as ambiguous);
  2. ChEMBL drug-mechanism action types, joined through the PDBeChem bioactivity report, used
     to fill gaps and to measure agreement with (1).
Partial agonists are kept as a separate descriptive class, not pooled with full agonists.

Unit of analysis: one volume per (PDB entry, pocket site) (median over that entry's pockets
assigned to the site), so multi-chain entries are not counted repeatedly. Redundant repeat
depositions were already collapsed by 04b_deduplicate_structures.py.

Needs: matched_pockets.csv, crypticity_index.csv, crypticity_index_with_family.csv,
ligand_descriptors.csv; network access (RCSB GraphQL, ChEMBL REST) only the first time --
labels are cached in data/external/nhr_entry_labels.csv (use --refresh to rebuild).

Output: results/nhr_agonist_antagonist_summary.txt
        results/nhr_site_comparison.csv
        results/nhr_entry_labels.csv
        results/figures/fig8_nhr_agonist_antagonist.png
"""
import argparse
import re
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy import stats

from common import load_config

MIN_PER_GROUP = 3
N_RESAMPLE = 500
SEED = 17
EPS = 1.0  # same epsilon as the crypticity index (A^3)
NAMES = {"P03372": "ERα", "P37231": "PPARγ", "P51449": "RORγt", "P04150": "GR", "P19793": "RXRα",
         "P62508": "ERRγ", "Q03181": "PPARδ", "Q07869": "PPARα", "O75469": "PXR"}
GROUP_OF = {"agonist": "A", "partial_agonist": "P", "antagonist": "N", "inverse_agonist": "N"}
CHEMBL_MAP = {"AGONIST": "agonist", "PARTIAL AGONIST": "partial_agonist", "ANTAGONIST": "antagonist",
              "INVERSE AGONIST": "inverse_agonist"}


# --------------------------------------------------------------------------- efficacy labels
def classify_text(text: str):
    s = text.lower()
    found = set()
    if "inverse agonist" in s:
        found.add("inverse_agonist")
    if "partial agonist" in s:
        found.add("partial_agonist")
    if re.search(r"antagonist|\bserm\b|\bserd\b|degrader", s):
        found.add("antagonist")
    rest = s.replace("inverse agonist", " ").replace("partial agonist", " ").replace("antagonist", " ")
    if "agonist" in rest:
        found.add("agonist")
    if not found:
        return None
    groups = {GROUP_OF[f] for f in found}
    if len(groups) > 1:
        return "ambiguous"
    return sorted(found)[0] if len(found) == 1 else ("inverse_agonist" if "inverse_agonist" in found else "antagonist")


def fetch_titles(pdb_ids):
    query = "query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title} struct_keywords{text pdbx_keywords}}}"
    rows = []
    for i in range(0, len(pdb_ids), 50):
        r = requests.post("https://data.rcsb.org/graphql", json={"query": query, "variables": {"ids": pdb_ids[i:i + 50]}}, timeout=60)
        r.raise_for_status()
        for e in r.json()["data"]["entries"] or []:
            if e:
                kw = e.get("struct_keywords") or {}
                rows.append({"pdb_id": e["rcsb_id"], "title": (e.get("struct") or {}).get("title") or "",
                             "keywords": f"{kw.get('text') or ''} {kw.get('pdbx_keywords') or ''}".strip()})
        time.sleep(0.2)
    return pd.DataFrame(rows)


def fetch_mechanisms(mol_ids):
    rows = []
    for i in range(0, len(mol_ids), 40):
        r = requests.get("https://www.ebi.ac.uk/chembl/api/data/mechanism.json",
                         params={"molecule_chembl_id__in": ",".join(mol_ids[i:i + 40]), "limit": 1000}, timeout=60)
        r.raise_for_status()
        rows += r.json().get("mechanisms", [])
        time.sleep(0.2)
    return pd.DataFrame(rows)


def build_labels(entries: pd.DataFrame, cfg) -> pd.DataFrame:
    titles = fetch_titles(sorted(entries["pdb_id"].unique()))
    lab = entries.merge(titles, on="pdb_id", how="left")
    lab["title"] = lab["title"].fillna("")
    lab["keywords"] = lab["keywords"].fillna("")
    lab["label_title"] = (lab["title"] + " " + lab["keywords"]).map(classify_text)
    lab["chembl_mol"], lab["label_chembl"] = None, None
    report = Path(cfg["bioactivity"]["report_path"]).with_name("bioactivity_report_simple.tsv")
    if report.exists():
        bio = pd.read_csv(report, sep="\t", usecols=["PDB", "LIG", "CHEMBL_ID_MOL"])
        bio["PDB"] = bio["PDB"].str.upper()
        lab["first_lig"] = lab["ligand_codes"].fillna("").str.split(";").str[0]
        lab = lab.drop(columns=["chembl_mol"]).merge(
            bio.drop_duplicates(["PDB", "LIG"]), left_on=["pdb_id", "first_lig"], right_on=["PDB", "LIG"], how="left"
        ).rename(columns={"CHEMBL_ID_MOL": "chembl_mol"}).drop(columns=["PDB", "LIG"])
        mols = sorted(lab["chembl_mol"].dropna().unique())
        mech = fetch_mechanisms(mols) if mols else pd.DataFrame()
        if len(mech):
            mech["lab"] = mech["action_type"].map(CHEMBL_MAP)
            mech = mech.dropna(subset=["lab"]).drop_duplicates("molecule_chembl_id")
            lab["label_chembl"] = lab["chembl_mol"].map(dict(zip(mech["molecule_chembl_id"], mech["lab"])))
    lab["label"] = lab["label_title"].where(lab["label_title"].isin(list(GROUP_OF)), lab["label_chembl"])
    lab["label_source"] = np.where(lab["label_title"].isin(list(GROUP_OF)), "title",
                                   np.where(lab["label_chembl"].notna(), "chembl", ""))
    lab["group"] = lab["label"].map(GROUP_OF)
    return lab


# --------------------------------------------------------------------------- statistics helpers
def rel_mad(x):
    x = np.asarray(x, float)
    m = np.median(x)
    return np.median(np.abs(x - m)) / m if m > 0 else np.nan


def equalised_dispersion(a, b, rng):
    k = min(len(a), len(b))
    da, db = [], []
    for _ in range(N_RESAMPLE):
        da.append(rel_mad(rng.choice(a, k, replace=False)))
        db.append(rel_mad(rng.choice(b, k, replace=False)))
    return np.nanmean(da), np.nanmean(db), k


def wilcoxon_line(label, d, expect=None):
    d = np.asarray(d, float)
    d = d[~np.isnan(d)]
    if len(d) < 6:
        return f"{label}: n={len(d)} sites -- too few to test."
    w, p = stats.wilcoxon(d)
    s = f"{label}: n={len(d)} sites, median difference={np.median(d):+.3f}, Wilcoxon signed-rank p={p:.3g}"
    if expect:
        s += f" ({'consistent with' if np.sign(np.median(d)) == expect[0] else 'opposite to'} expectation: {expect[1]})"
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true", help="re-fetch efficacy labels (RCSB + ChEMBL)")
    ap.add_argument("--family", default="nuclear hormone receptor", help="substring of the UniProt family string")
    args = ap.parse_args()

    cfg = load_config()
    res, root = Path(cfg["paths"]["results"]), Path(cfg["paths"]["results"]).parent
    try:
        matched = pd.read_csv(f"{cfg['paths']['pockets']}/matched_pockets.csv",
                              usecols=["pdb_id", "uniprot", "is_apo", "ligand_codes", "resolution", "pocket_site_id", "volume"])
        crypt = pd.read_csv(res / "crypticity_index_with_family.csv")
        lig = pd.read_csv(f"{cfg['paths']['ligands']}/ligand_descriptors.csv", usecols=["het_code", "heavy_atoms"])
    except FileNotFoundError as e:
        sys.exit(f"Missing input ({e}). Run the pipeline through 07/10 first.")

    fam = crypt["protein_family"].fillna("").str.contains(args.family, case=False)
    prots = sorted(crypt.loc[fam, "uniprot"].unique())
    holo = matched[(~matched["is_apo"]) & matched["uniprot"].isin(prots)].copy()
    holo["pdb_id"] = holo["pdb_id"].str.upper()
    entries = holo[["pdb_id", "uniprot", "ligand_codes"]].drop_duplicates("pdb_id")

    cache = root / "data/external/nhr_entry_labels.csv"
    if cache.exists() and not args.refresh:
        labels = pd.read_csv(cache, keep_default_na=False, na_values=[""])
    else:
        labels = build_labels(entries, cfg)
        cache.parent.mkdir(parents=True, exist_ok=True)
        labels.to_csv(cache, index=False)
    labels.to_csv(res / "nhr_entry_labels.csv", index=False)

    L = ["=" * 72, f"Hypothesis F: ligand efficacy vs. pocket plasticity ({args.family})", "=" * 72, "",
         f"{len(prots)} proteins, {len(entries)} holo entries in the deduplicated dataset."]
    n_lab = labels["group"].notna().sum()
    L.append(f"Efficacy label available for {n_lab} entries ({100 * n_lab / len(labels):.0f}%): "
             f"{(labels['label_source'] == 'title').sum()} from the PDB title/keywords, "
             f"{(labels['label_source'] == 'chembl').sum()} added from ChEMBL mechanisms.")
    both = labels[labels["label_title"].isin(list(GROUP_OF)) & labels["label_chembl"].isin(list(GROUP_OF))]
    if len(both):
        agree = (both["label_title"].map(GROUP_OF) == both["label_chembl"].map(GROUP_OF)).mean()
        L.append(f"Title vs. ChEMBL labels agree (A/P/N class) for {100 * agree:.0f}% of the {len(both)} entries having both.")
    L += ["", "Labelled entries per protein and class (A=agonist, P=partial agonist, N=antagonist/inverse agonist):"]
    ct = labels.dropna(subset=["group"]).groupby(["uniprot", "group"]).size().unstack(fill_value=0)
    ct.index = [f"{NAMES.get(u, u)} ({u})" for u in ct.index]
    L.append(ct.to_string())

    # ---- observation table: one volume per (entry, site)
    obs = holo.merge(labels[["pdb_id", "group"]], on="pdb_id").dropna(subset=["group"])
    obs = obs.groupby(["pdb_id", "uniprot", "pocket_site_id", "group"], as_index=False).agg(
        volume=("volume", "median"), resolution=("resolution", "first"), ligand_codes=("ligand_codes", "first"))
    ha = dict(zip(lig["het_code"], lig["heavy_atoms"]))
    obs["heavy_atoms"] = obs["ligand_codes"].fillna("").map(
        lambda s: max([ha.get(c, np.nan) for c in s.split(";") if c] or [np.nan]))
    obs = obs.merge(crypt[["pocket_site_id", "v_apo_min"]], on="pocket_site_id")

    # ---- site-level comparison
    rng = np.random.default_rng(SEED)
    rows = []
    for site, g in obs.groupby("pocket_site_id"):
        a, n = g.loc[g["group"] == "A", "volume"].values, g.loc[g["group"] == "N", "volume"].values
        if len(a) < MIN_PER_GROUP or len(n) < MIN_PER_GROUP:
            continue
        da, dn, k = equalised_dispersion(a, n, rng)
        va = g["v_apo_min"].iloc[0]
        rows.append({"pocket_site_id": site, "uniprot": g["uniprot"].iloc[0], "n_agonist": len(a), "n_antagonist": len(n),
                     "n_partial": int((g["group"] == "P").sum()), "median_vol_agonist": np.median(a),
                     "median_vol_antagonist": np.median(n), "rel_mad_agonist_eq": da, "rel_mad_antagonist_eq": dn,
                     "resample_k": k, "ci_agonist": np.median(a) / (va + EPS), "ci_antagonist": np.median(n) / (va + EPS),
                     "v_apo_min": va, "mw_p_volume": stats.mannwhitneyu(a, n).pvalue})
    sites = pd.DataFrame(rows)
    sites.to_csv(res / "nhr_site_comparison.csv", index=False)
    L += ["", f"Sites populated by >= {MIN_PER_GROUP} agonist AND >= {MIN_PER_GROUP} antagonist-like entries: {len(sites)} "
          f"(proteins: {', '.join(NAMES.get(u, u) for u in sorted(sites['uniprot'].unique())) if len(sites) else 'none'})."]

    fig_ok = len(sites) >= 1
    if len(sites):
        L += ["", "--- Site-level paired tests (agonist vs. antagonist-like) ---",
              wilcoxon_line("F1 volume   (antagonist - agonist, A^3)", sites["median_vol_antagonist"] - sites["median_vol_agonist"]),
              wilcoxon_line("F2 plasticity (relative MAD, agonist - antagonist; expected < 0)",
                            sites["rel_mad_agonist_eq"] - sites["rel_mad_antagonist_eq"], (-1, "agonist-bound pockets vary less")),
              wilcoxon_line("F3 crypticity (CI_antagonist - CI_agonist; expected > 0)",
                            sites["ci_antagonist"] - sites["ci_agonist"], (1, "antagonist-like ligands open the pocket more"))]
        per_prot = sites.groupby("uniprot").agg(sites=("pocket_site_id", "size"),
            vol_A=("median_vol_agonist", "median"), vol_N=("median_vol_antagonist", "median"),
            relMAD_A=("rel_mad_agonist_eq", "median"), relMAD_N=("rel_mad_antagonist_eq", "median"),
            CI_A=("ci_agonist", "median"), CI_N=("ci_antagonist", "median")).round(3)

        def _p(g, a, b):
            d = (g[a] - g[b]).dropna()
            return round(stats.wilcoxon(d).pvalue, 4) if len(d) >= 6 and (d != 0).any() else np.nan
        grp = sites.groupby("uniprot")
        per_prot["p_vol"] = grp.apply(lambda g: _p(g, "median_vol_antagonist", "median_vol_agonist"))
        per_prot["p_relMAD"] = grp.apply(lambda g: _p(g, "rel_mad_agonist_eq", "rel_mad_antagonist_eq"))
        per_prot["p_CI"] = grp.apply(lambda g: _p(g, "ci_antagonist", "ci_agonist"))
        per_prot.index = [f"{NAMES.get(u, u)} ({u})" for u in per_prot.index]
        L += ["", "Per-protein medians over qualifying sites (p = within-protein Wilcoxon over sites; sites of one",
              "protein are not independent, so read these as descriptive):", per_prot.to_string()]

        # ---- orthosteric site per protein: entry-level view
        L += ["", "--- Most populated site per protein (entry-level Mann-Whitney, A vs N volume) ---"]
        top = obs[obs["group"].isin(["A", "P", "N"])].groupby(["uniprot", "pocket_site_id"]).size().reset_index(name="n")
        top = top.sort_values("n", ascending=False).drop_duplicates("uniprot")
        top_sites = {}
        for _, r in top.iterrows():
            g = obs[obs["pocket_site_id"] == r["pocket_site_id"]]
            a, p_, n = (g.loc[g["group"] == x, "volume"].values for x in "APN")
            top_sites[r["uniprot"]] = r["pocket_site_id"]
            if len(a) >= MIN_PER_GROUP and len(n) >= MIN_PER_GROUP:
                L.append(f"{NAMES.get(r['uniprot'], r['uniprot']):6s} {r['pocket_site_id']}: median volume A={np.median(a):.0f} (n={len(a)}), "
                         f"P={np.median(p_) if len(p_) else float('nan'):.0f} (n={len(p_)}), N={np.median(n):.0f} (n={len(n)}) A^3; "
                         f"MWU p={stats.mannwhitneyu(a, n).pvalue:.3g}")

        # ---- mixed model on entry x site observations
        try:
            import statsmodels.formula.api as smf
            d = obs[obs["group"].isin(["A", "N"]) & obs["pocket_site_id"].isin(sites["pocket_site_id"])].copy()
            d["logvol"], d["N"] = np.log(d["volume"]), (d["group"] == "N").astype(int)
            m1 = smf.mixedlm("logvol ~ N", d, groups=d["pocket_site_id"]).fit(reml=True)
            L += ["", "--- Mixed model, log volume ~ class + (1 | site) ---",
                  f"antagonist-like vs agonist: coefficient={m1.params['N']:+.3f} log units "
                  f"(x{np.exp(m1.params['N']):.2f}), p={m1.pvalues['N']:.3g}, {len(d)} observations / {d['pocket_site_id'].nunique()} sites"]
            d2 = d.dropna(subset=["heavy_atoms", "resolution"]).copy()
            for c in ("heavy_atoms", "resolution"):
                d2[c + "_z"] = (d2[c] - d2[c].mean()) / d2[c].std()
            m2 = smf.mixedlm("logvol ~ N + heavy_atoms_z + resolution_z", d2, groups=d2["pocket_site_id"]).fit(reml=True)
            L.append(f"adjusted for ligand heavy-atom count and resolution: class coefficient={m2.params['N']:+.3f} "
                     f"(x{np.exp(m2.params['N']):.2f}), p={m2.pvalues['N']:.3g}; heavy-atom coefficient="
                     f"{m2.params['heavy_atoms_z']:+.3f}, p={m2.pvalues['heavy_atoms_z']:.3g} (n={len(d2)})")
        except Exception as e:
            L.append(f"[mixed model skipped: {e}]")

        # ---- partial-agonist sensitivity (agonist-like = A + P vs N)
        al = obs[obs["group"].isin(["A", "P", "N"])].copy()
        al["g2"] = np.where(al["group"] == "N", "N", "A")
        rows2 = []
        for site, g in al.groupby("pocket_site_id"):
            a, n = g.loc[g["g2"] == "A", "volume"].values, g.loc[g["g2"] == "N", "volume"].values
            if len(a) >= MIN_PER_GROUP and len(n) >= MIN_PER_GROUP:
                rows2.append((np.median(n) - np.median(a)))
        L += ["", "--- Sensitivity: partial agonists pooled with agonists ---",
              wilcoxon_line("F1 volume (antagonist - [agonist+partial], A^3)", rows2)]

    L += ["", "Reading guide: only proteins with enough structures of BOTH classes can be tested; labels come from PDB",
          "titles (not pharmacology), ligand size is confounded with efficacy class, and sites within one protein are",
          "not independent -- see the mixed model and the per-protein table before drawing family-wide conclusions."]
    report = "\n".join(L)
    print(report)
    (res / "nhr_agonist_antagonist_summary.txt").write_text(report, encoding="utf-8")

    # ---- figure
    if fig_ok:
        plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "font.size": 9})
        fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.9))
        col = {"A": "#2A9D8F", "P": "#E9C46A", "N": "#C44E52"}
        pos, ticks, labs = 0, [], []
        for u, s in top_sites.items():
            g = obs[obs["pocket_site_id"] == s]
            if (g["group"] == "A").sum() < MIN_PER_GROUP or (g["group"] == "N").sum() < MIN_PER_GROUP:
                continue
            for x in "APN":
                v = g.loc[g["group"] == x, "volume"].values
                if len(v):
                    ax[0].boxplot(v, positions=[pos], widths=0.7, showfliers=False, patch_artist=True,
                                  boxprops=dict(facecolor=col[x], alpha=0.6), medianprops=dict(color="k"))
                    ax[0].scatter(pos + rng.uniform(-0.2, 0.2, len(v)), v, s=6, color="k", alpha=0.35, zorder=3)
                pos += 1
            ticks.append(pos - 2)
            labs.append(NAMES.get(u, u))
            pos += 1
        ax[0].set_xticks(ticks); ax[0].set_xticklabels(labs)
        ax[0].set_ylabel("pocket volume (Å$^3$)")
        ax[0].set_title("Most populated site per protein", fontsize=9, fontweight="bold")
        for x, nme in (("A", "agonist"), ("P", "partial"), ("N", "antagonist-like")):
            ax[0].bar(0, 0, color=col[x], alpha=0.6, label=nme)
        ax[0].legend(frameon=False, fontsize=7)
        pc = {u: c for u, c in zip(sorted(sites["uniprot"].unique()), ["#4C72B0", "#DD8452", "#55A868", "#8172B2", "#937860"])}
        for u, g in sites.groupby("uniprot"):
            ax[1].scatter(g["median_vol_agonist"], g["median_vol_antagonist"], s=14, color=pc[u], alpha=0.8, label=NAMES.get(u, u))
            ax[2].scatter(g["rel_mad_agonist_eq"], g["rel_mad_antagonist_eq"], s=14, color=pc[u], alpha=0.8)
        for a_, ttl, xl, yl in ((ax[1], "Median volume per site", "agonist-bound (Å$^3$)", "antagonist-like (Å$^3$)"),
                                (ax[2], "Plasticity per site (relative MAD)", "agonist-bound", "antagonist-like")):
            lo = min(a_.get_xlim()[0], a_.get_ylim()[0]); hi = max(a_.get_xlim()[1], a_.get_ylim()[1])
            a_.plot([lo, hi], [lo, hi], color="grey", lw=0.8, ls="--"); a_.set_xlim(lo, hi); a_.set_ylim(lo, hi)
            a_.set_title(ttl, fontsize=9, fontweight="bold"); a_.set_xlabel(xl); a_.set_ylabel(yl)
        ax[1].legend(frameon=False, fontsize=7)
        fig.tight_layout()
        out = res / "figures" / "fig8_nhr_agonist_antagonist.png"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
