# Cryptic Pocket Prediction Pipeline

Tests the hypothesis: structural variability across independently solved
crystal structures of the same human protein predicts cryptic-pocket
formation, and that crypticity in turn predicts ligand presence, size,
chemical diversity, and affinity.

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
12_make_figures.py             5 figures documenting the hypothesis results ->
                                 results/figures/*.png
13_join_bioactivity.py         Hypothesis D (affinity): joins the max-volume-
                                 holo ligand against real ChEMBL bioactivity
                                 data (external PDBeChem report, see config)
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
   an external PDBeChem/ChEMBL bioactivity report (18.7% site coverage).
   Result so far: rho=-0.038, p=0.018 on n=3,908 — statistically significant
   but a negligible effect size, and even the sign runs the "wrong" way
   (very weakly, cryptic sites trend toward *less* potent ligands). Don't
   read this as a real relationship; see `results/affinity_correlation.txt`.
4. **Chemical diversity (hypothesis C) needs many holo structures per site**
   to have any power — most pocket sites will only ever show 1 bound ligand
   across the PDB. This arm is realistically only testable after a broad
   (many-protein) run, not a handful of hand-picked targets.
5. **fpocket cavity detection is somewhat parameter-sensitive** — its default
   alpha-sphere clustering can merge or split what a human would call one
   pocket. Cross-validating against a second, independent pocket detector
   (DoGSite3 was tried and dropped — see below) would help sanity-check the
   crypticity calls, but nothing is currently wired up for this.
6. **Entries only available as mmCIF with multi-letter chain IDs (e.g. "AAA")
   are skipped** in `03_download_and_clean.py` — legacy PDB format caps chain
   IDs at 1 character, so `PDBIO.save` throws. Affected a handful of the 53
   CA1 entries in the smoke test. Fixable by writing cleaned structures out
   as mmCIF instead of PDB if this starts discarding too much data on a
   broader run.
