# Cryptic Pocket Prediction Pipeline

Tests the hypothesis: structural variability across independently solved
crystal structures of the same human protein predicts cryptic-pocket
formation, and that crypticity in turn predicts ligand presence, size,
chemical diversity, and affinity — plus (hypothesis E) that a standard
pocket detector's own apo-state druggability score underrates sites that
later turn out to be cryptic, the premise that motivated this whole project.

    apo (small pocket) --------> holo (large pocket, ligand bound)
    V_apo                        V_liganded
    crypticity_index = V_liganded / (V_apo + epsilon)
    delta_V_pocket   = V_liganded - V_apo

Pipeline stages (run in order from `scripts/`):

```
01_fetch_pdb_index.py          RCSB Search API, grouped by UniProt accession ->
                                 all human/X-ray/<=2.5A entries per protein
02_classify_apo_holo.py        Apo vs holo classification (blacklist-filtered
                                 ligand HET codes) -> keep proteins with both
03_download_and_clean.py       Download PDB files, strip to protein-only atoms
04_run_fpocket.py               fpocket on every structure -> every pocket's
                                 volume/druggability/centroid/lining residues
04b_deduplicate_structures.py  Redundancy reduction: collapses repeat holo
                                 depositions of the same protein+ligand complex
                                 to one best-resolution representative (apo
                                 structures untouched — their variability IS
                                 the signal). 05+ prefer this output when present.
05_align_and_match_pockets.py  Bio.PDB CA superposition onto a reference
                                 structure per protein; greedy-cluster pockets
                                 across structures into pocket_site_id's
06_compute_crypticity.py       Per-site delta_V_pocket + crypticity_index
07_fetch_ligand_descriptors.py RCSB chem-comp dictionary + RDKit descriptors
                                 (heavy atoms, MW, logP, TPSA, Murcko scaffold)
08_correlate_crypticity.py     Ligand presence (Mann-Whitney U + logistic
                                 regression), size and chemical-diversity
                                 correlations (Spearman) against crypticity_index
10_fetch_protein_families.py   UniProt family/keyword join, for stratifying
                                 results by protein class (kinases vs. not, etc.)
11_stratify_by_family.py       Per-family crypticity stats + within-family
                                 diversity correlation; overall Kruskal-Wallis
                                 across families (formalizes the ad hoc check)
12_make_figures.py             5-7 figures documenting the hypothesis results ->
                                 results/figures/*.png (fig6/fig7 need 13/14
                                 to have been run first; skipped otherwise)
13_join_bioactivity.py         Hypothesis D (affinity): joins the max-volume-
                                 holo ligand against real ChEMBL bioactivity
                                 data (external PDBeChem report, see config)
14_druggability_analysis.py    Hypothesis E: fpocket's own apo-state
                                 druggability score (already computed in 06,
                                 never analyzed) vs. crypticity_index
```

## Setup

```bash
pip install -r requirements.txt
```

External tools used (already present on this machine, per `common.py`):
- `fpocket` — resolved on PATH, else run through WSL (confirmed at
  `/usr/local/bin/fpocket` in this machine's WSL). Same wrapper as
  `NEW_IDEA/hotspot_mimicry_pipeline`.
- `rdkit`, `biopython` — via the base conda/pip environment.
- OpenBabel and PyMOL are available on this machine (`Documents/pymol-open-source`,
  native `obabel.exe`) but are **not currently used** — pocket volumes come
  straight from fpocket, and structure superposition uses `Bio.PDB.Superimposer`
  rather than PyMOL's `align`/`cealign`. See "Limitations" below for where
  swapping in PyMOL's structural alignment would help.

## Configuration (`config.yaml`)

- `rcsb.target_uniprot_ids`: leave `[]` to run the full discovery crawl over
  all matching human X-ray structures (capped at `rcsb.max_groups` UniProt
  groups, ranked by structure count), or fill in specific UniProt accessions
  to fast-path a hand-picked target list (recommended for a first run — the
  full crawl currently matches ~5,900 UniProt groups / ~60,000 entries).
- `selection`: minimum apo/holo counts required to keep a protein.
  `solvent_blacklist`: HET codes never counted as a real ligand.
- `pocket_matching`: thresholds for calling two structures' pockets "the
  same site" after superposition.

## Known limitations (be aware of these before trusting the numbers)

1. **Cross-structure alignment matches by PDB residue number**, not by
   sequence alignment. This is fine for most well-behaved entries of a
   protein solved by the same group/numbering convention, but breaks for
   entries with inconsistent numbering, engineered constructs, or multiple
   isoforms — those structures silently fall back to `aligned=False` in
   `matched_pockets.csv` (their pocket keeps its raw, unaligned centroid, so
   it likely won't cluster into any real site). Upgrading to a proper
   sequence alignment (or PyMOL's `cealign`, which is structure- rather than
   numbering-based) would fix this and is the single highest-value follow-up.
2. **Pocket clustering is a greedy single-pass heuristic**, not a proper
   graph/hierarchical clustering — cluster order depends on input order.
   Fine for scaffolding, worth revisiting if site assignment looks noisy.
3. **Affinity (hypothesis D)**: computed in `13_join_bioactivity.py` against
   an external PDBeChem/ChEMBL bioactivity report (19.3% site coverage,
   n=3,654). Result: rho=-0.031, p=0.060 — not even nominally significant
   after redundancy reduction (see item 7 below; it was rho=-0.038, p=0.018
   *before* deduplication, so the earlier borderline "significance" looks
   like it was partly a pseudo-replication artifact). Read this as a null
   result, not a real relationship; see `results/affinity_correlation.txt`.
4. **Chemical diversity (hypothesis C, `n_distinct_scaffolds`) still
   conflates two things**: how chemically accommodating a site genuinely is,
   and how many independent structure-based drug discovery campaigns have
   targeted that protein. Redundancy reduction (item 7) removes duplicate
   depositions of the *same* complex, but not the effect of a protein having
   genuinely been co-crystallized with many *different* real ligands across
   many separate studies — beta-2-microglobulin (877 sites even after
   dedup, reflecting hundreds of distinct, non-redundant peptide-MHC
   complexes) is the clearest example. Normalizing `n_distinct_scaffolds` by
   the number of independent holo structures at a site (or by a research-
   intensity proxy) would be needed to properly separate "genuinely
   promiscuous pocket" from "heavily studied protein."
5. **fpocket cavity detection is somewhat parameter-sensitive** — its default
   alpha-sphere clustering can merge or split what a human would call one
   pocket. Cross-validating against a second, independent pocket detector
   (DoGSite3 was tried and dropped — see below) would help sanity-check the
   crypticity calls, but nothing is currently wired up for this.
6. **Entries only available as mmCIF with multi-letter chain IDs (e.g. "AAA")
   are skipped** in `03_download_and_clean.py` — legacy PDB format caps chain
   IDs at 1 character, so `PDBIO.save` throws. A minority of entries across
   the full run hit this. Fixable by writing cleaned structures out as
   mmCIF instead of PDB if this starts discarding too much data on a
   broader run.
7. **Redundancy reduction (`04b_deduplicate_structures.py`) is applied but
   partial.** It collapses exact repeat depositions of the same protein
   sharing an identical bound-ligand set to one best-resolution
   representative — this removed 21.3% of holo structures (4,827 of 22,661)
   and roughly halved the most extreme outlier (beta-2-microglobulin: 1,387
   → 877 pocket sites). It does **not** address: (a) non-independence
   between *distinct* sites of the same protein — a protein contributing
   dozens of separate `pocket_site_id`s still inflates the pooled sample
   size relative to the number of truly independent proteins; or (b) the
   research-intensity confound in item 4 above. Pooled statistical tests
   (Mann-Whitney U, Spearman, Kruskal-Wallis throughout `08`/`11`) still
   treat every site as an independent observation, which isn't strictly
   true. A mixed-effects model with protein as a random effect is the
   proper fix and hasn't been done. That said, the fact that this step
   pushed hypothesis D from nominally significant to non-significant (item
   3) while leaving hypotheses A and C essentially unchanged suggests it's
   doing real, not cosmetic, work — not just an unverified caveat.
8. **Druggability score (hypothesis E) is heavily zero-inflated and not
   fully independent of crypticity_index.** `druggability_apo_min_v` /
   `druggability_holo_max_v` are fpocket's own score for the specific
   matched pocket at a site, not necessarily a structure's single
   top-ranked pocket, so most sites — cryptic or not — score near 0
   (dataset median 0.001). The CI≥5 vs. CI<2 comparison is real
   (p=4.2×10⁻⁴, mean 0.015 vs. 0.036) but is population-level evidence
   that the phenomenon exists, not a usable per-site predictor. It's also
   not fully independent of Hypothesis A/crypticity_index itself, since
   fpocket's druggability score is partly a function of the same pocket
   volume.

## Data to be downloaded into /data folder

Use the fowllowing Zenodo [link](https://doi.org/10.5281/zenodo.23037984)
