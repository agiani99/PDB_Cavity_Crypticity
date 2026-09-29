"""
Discover UniProt-grouped sets of human, X-ray, resolution <= cutoff PDB
entries. Either restricts to `config.yaml`'s rcsb.target_uniprot_ids (fast
path for a hand-picked target list) or crawls the full RCSB grouping and
keeps the `max_groups` groups with the most member structures (broad
discovery, capped for tractability).

Output: data/index/pdb_groups.csv
    uniprot, pdb_id, entity_id, resolution
"""
import sys
import pandas as pd
from common import load_config, rcsb_group_by_uniprot, rcsb_entry_details


def main():
    cfg = load_config()
    rc = cfg["rcsb"]

    print("Querying RCSB Search API (human, X-ray, resolution <= "
          f"{rc['resolution_max']} A), grouped by UniProt accession...")
    groups = rcsb_group_by_uniprot(cfg)
    print(f"Fetched {len(groups)} UniProt groups.")

    target_ids = set(rc.get("target_uniprot_ids") or [])
    if target_ids:
        groups = [g for g in groups if g["identifier"] in target_ids]
        print(f"Restricted to {len(groups)} target UniProt accessions.")
        missing = target_ids - {g["identifier"] for g in groups}
        if missing:
            print(f"[warn] no matching structures found for: {sorted(missing)}")
    else:
        groups = [g for g in groups if g["count"] >= cfg["selection"]["min_structures_per_protein"]]
        groups.sort(key=lambda g: g["count"], reverse=True)
        groups = groups[: rc["max_groups"]]
        print(f"Discovery mode: kept top {len(groups)} groups by structure count "
              f"(min {cfg['selection']['min_structures_per_protein']} structures each).")

    rows = []
    for g in groups:
        uniprot = g["identifier"]
        for member in g["result_set"]:
            pdb_id, _, entity_id = member["identifier"].partition("_")
            rows.append({"uniprot": uniprot, "pdb_id": pdb_id, "entity_id": entity_id})

    if not rows:
        sys.exit("No groups/members found — check config.yaml [rcsb] settings.")

    df = pd.DataFrame(rows).drop_duplicates(subset=["uniprot", "pdb_id"])
    print(f"Resolving resolution + ligand HET codes for {df['pdb_id'].nunique()} entries...")
    details = rcsb_entry_details(df["pdb_id"].tolist(), cfg)
    df["resolution"] = df["pdb_id"].map(lambda e: details.get(e, {}).get("resolution"))
    df["het_codes"] = df["pdb_id"].map(lambda e: ";".join(details.get(e, {}).get("het_codes", [])))

    out_path = f"{cfg['paths']['index']}/pdb_groups.csv"
    df.to_csv(out_path, index=False)
    print(f"Wrote {len(df)} entries across {df['uniprot'].nunique()} proteins -> {out_path}")


if __name__ == "__main__":
    main()
