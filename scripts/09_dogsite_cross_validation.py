"""
Cross-validate fpocket-derived crypticity calls against DoGSite3
(proteins.plus) on a small subset — the reference apo structure and the
max-volume holo structure for each of the top-N most cryptic sites.

STATUS: blocked. Every payload shape tried live against
https://proteins.plus/api/dogsite3_rest returned
'400 Invalid number of parameters or incorrect parameter name', including
the JSON/form, flat/nested-under-"dogsite3" variants documented in various
proteins.plus examples. Register at https://proteins.plus/api for an API
key — the confirmation email includes a working curl example — then fix
`run_dogsite3()` in common.py accordingly. This script is otherwise wired
up to run once that function works.
"""
import sys
import pandas as pd
from common import load_config, run_dogsite3


def main():
    cfg = load_config()
    try:
        crypticity = pd.read_csv(f"{cfg['paths']['results']}/crypticity_index.csv")
    except FileNotFoundError:
        sys.exit("Missing results/crypticity_index.csv. Run 06_compute_crypticity.py first.")

    top_n = 10
    subset = crypticity.head(top_n)
    rows = []
    for _, row in subset.iterrows():
        for pdb_id, label in [(row["holo_pdb_at_max"], "holo_max")]:
            try:
                result = run_dogsite3(pdb_id, cfg)
                rows.append({"pocket_site_id": row["pocket_site_id"], "pdb_id": pdb_id,
                             "label": label, **result})
            except NotImplementedError as e:
                sys.exit(str(e))
            except Exception as e:
                print(f"[skip] DoGSite3 failed for {pdb_id}: {e}")

    out = pd.DataFrame(rows)
    out_path = f"{cfg['paths']['results']}/dogsite3_cross_validation.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} DoGSite3 results -> {out_path}")


if __name__ == "__main__":
    main()
