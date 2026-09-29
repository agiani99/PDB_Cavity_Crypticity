"""
Join UniProt protein-family membership + broad functional keywords onto the
per-site crypticity results, so downstream analysis can stratify by protein
class (e.g. "are cryptic pockets more common in kinases than hydrolases?")
without a join at analysis time.

Output: results/crypticity_index_with_family.csv
    ...all crypticity_index.csv columns..., protein_family, keywords
"""
import sys
import pandas as pd
from common import load_config, uniprot_protein_family


def main():
    cfg = load_config()
    in_path = f"{cfg['paths']['results']}/crypticity_index.csv"
    try:
        df = pd.read_csv(in_path)
    except FileNotFoundError:
        sys.exit(f"Missing {in_path}. Run 06_compute_crypticity.py first.")

    uniprot_ids = df["uniprot"].dropna().unique().tolist()
    print(f"Fetching UniProt family/keyword data for {len(uniprot_ids)} accessions...")
    families = uniprot_protein_family(uniprot_ids, cfg)

    df["protein_family"] = df["uniprot"].map(lambda u: families.get(u, {}).get("protein_family"))
    df["keywords"] = df["uniprot"].map(lambda u: families.get(u, {}).get("keywords"))

    missing = [u for u in uniprot_ids if u not in families]
    if missing:
        print(f"[warn] no UniProt record found for: {missing}")

    out_path = f"{cfg['paths']['results']}/crypticity_index_with_family.csv"
    df.to_csv(out_path, index=False)
    print(f"Wrote {len(df)} rows -> {out_path}")

    with_family = df["protein_family"].notna().sum()
    print(f"{with_family}/{len(df)} sites resolved to a protein family, "
          f"covering {df.loc[df['protein_family'].notna(), 'uniprot'].nunique()} proteins.")


if __name__ == "__main__":
    main()
