"""Shared helpers for the cryptic-pocket pipeline."""
import shutil
import subprocess
import time
from pathlib import Path

import requests
import yaml

ROOT = Path(__file__).resolve().parent.parent


def load_config(path: str = None) -> dict:
    cfg_path = Path(path) if path else ROOT / "config.yaml"
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)
    for key, val in cfg["paths"].items():
        cfg["paths"][key] = str(ROOT / val)
        Path(cfg["paths"][key]).mkdir(parents=True, exist_ok=True)
    return cfg


# ---------------------------------------------------------------------------
# WSL / fpocket plumbing (same approach as hotspot_mimicry_pipeline/common.py)
# ---------------------------------------------------------------------------

def win_to_wsl_path(win_path: str) -> str:
    """Convert a native Windows path to its /mnt/c/... WSL equivalent."""
    normalized = str(win_path).replace("\\", "/")
    result = subprocess.run(
        ["wsl.exe", "wslpath", "-a", normalized],
        check=True, capture_output=True, text=True,
    )
    return result.stdout.strip()


def run_wsl(binary_win_path: str, args: list):
    """Run a Linux-only binary (native Windows path) inside WSL, converting
    every path-shaped argument to /mnt/c/... form first."""
    wsl_binary = win_to_wsl_path(binary_win_path) if ("\\" in binary_win_path or "/" in binary_win_path) else binary_win_path
    wsl_args = [
        win_to_wsl_path(a) if isinstance(a, str) and ("\\" in a or "/" in a) else str(a)
        for a in args
    ]
    return subprocess.run(
        ["wsl.exe", "-e", wsl_binary, *wsl_args],
        check=True, capture_output=True, text=True,
    )


def run_fpocket(binary: str, pdb_path: str):
    """Run fpocket on `pdb_path`. Uses the native binary if on PATH,
    otherwise shells into WSL (fpocket is Linux-only)."""
    if shutil.which(binary):
        subprocess.run([binary, "-f", pdb_path], check=True, capture_output=True)
        return
    wsl_path = win_to_wsl_path(pdb_path)
    subprocess.run(
        ["wsl.exe", "-e", binary, "-f", wsl_path],
        check=True, capture_output=True,
    )


def fpocket_output_dir(pdb_path) -> Path:
    return Path(f"{Path(pdb_path).with_suffix('')}_out")


# ---------------------------------------------------------------------------
# RCSB Search API + Data API (confirmed working against the live endpoints)
# ---------------------------------------------------------------------------

def rcsb_group_by_uniprot(cfg: dict) -> list[dict]:
    """Query the RCSB Search API for human X-ray structures at or below the
    configured resolution cutoff, grouped by matching UniProt accession.
    Returns a list of {identifier, count, result_set:[{identifier: "PDBID_entity"}...]}.
    """
    rc = cfg["rcsb"]
    query = {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "exptl.method", "operator": "exact_match", "value": rc["method"]}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less_or_equal", "value": rc["resolution_max"]}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.id",
                "operator": "exact_match", "value": str(rc["taxonomy_id"])}},
        ],
    }
    page = rc["group_page_size"]
    groups, start = [], 0
    while True:
        payload = {
            "query": query,
            "return_type": "polymer_entity",
            "request_options": {
                "group_by": {"aggregation_method": "matching_uniprot_accession"},
                "group_by_return_type": "groups",
                "paginate": {"start": start, "rows": page},
            },
        }
        r = requests.post(rc["search_url"], json=payload, timeout=30)
        r.raise_for_status()
        data = r.json()
        chunk = data.get("group_set", [])
        groups.extend(chunk)
        total = data.get("group_by_count", len(groups))
        start += page
        time.sleep(rc["request_pause_s"])
        if start >= total or not chunk:
            break
    return groups


def rcsb_entry_details(entry_ids: list[str], cfg: dict) -> dict:
    """Batch-fetch resolution + non-polymer (ligand) HET codes for entry IDs
    via the RCSB Data API GraphQL endpoint. Returns {entry_id: {resolution, het_codes}}."""
    rc = cfg["rcsb"]
    query = """
    query($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        rcsb_entry_info { resolution_combined }
        nonpolymer_entities { pdbx_entity_nonpoly { comp_id name } }
      }
    }
    """
    out = {}
    batch_size = rc["entry_batch_size"]
    unique_ids = sorted(set(entry_ids))
    for i in range(0, len(unique_ids), batch_size):
        batch = unique_ids[i:i + batch_size]
        r = requests.post(rc["data_graphql_url"],
                           json={"query": query, "variables": {"ids": batch}},
                           timeout=30)
        r.raise_for_status()
        entries = r.json().get("data", {}).get("entries") or []
        for e in entries:
            if e is None:
                continue
            res = e.get("rcsb_entry_info", {}).get("resolution_combined") or [None]
            hets = [
                (ent.get("pdbx_entity_nonpoly") or {}).get("comp_id")
                for ent in (e.get("nonpolymer_entities") or [])
            ]
            out[e["rcsb_id"]] = {
                "resolution": res[0],
                "het_codes": [h for h in hets if h],
            }
        time.sleep(rc["request_pause_s"])
    return out


def rcsb_chem_comp_descriptors(comp_ids: list[str], cfg: dict) -> dict:
    """Batch-fetch formula / weight / SMILES for HET (chemical component) codes."""
    rc = cfg["rcsb"]
    query = """
    query($ids: [String!]!) {
      chem_comps(comp_ids: $ids) {
        rcsb_id
        chem_comp { formula formula_weight name }
        rcsb_chem_comp_descriptor { SMILES }
      }
    }
    """
    out = {}
    batch_size = rc["entry_batch_size"]
    unique_ids = sorted(set(comp_ids))
    for i in range(0, len(unique_ids), batch_size):
        batch = unique_ids[i:i + batch_size]
        r = requests.post(rc["data_graphql_url"],
                           json={"query": query, "variables": {"ids": batch}},
                           timeout=30)
        r.raise_for_status()
        comps = r.json().get("data", {}).get("chem_comps") or []
        for c in comps:
            if c is None:
                continue
            cc = c.get("chem_comp") or {}
            desc = c.get("rcsb_chem_comp_descriptor") or {}
            out[c["rcsb_id"]] = {
                "formula": cc.get("formula"),
                "formula_weight": cc.get("formula_weight"),
                "name": cc.get("name"),
                "smiles": desc.get("SMILES"),
            }
        time.sleep(rc["request_pause_s"])
    return out


def download_pdb_file(entry_id: str, dest_dir: str, cfg: dict) -> Path:
    """Download an entry as legacy .pdb, falling back to .cif for entries
    RCSB only serves in mmCIF format (newer/large depositions). Returns the
    path with whichever suffix was actually fetched; callers must be able
    to parse either (Bio.PDB's PDBParser/MMCIFParser both work)."""
    for ext in ("pdb", "cif"):
        dest = Path(dest_dir) / f"{entry_id}.{ext}"
        if dest.exists():
            return dest
        url = f"{cfg['rcsb']['files_base_url']}/{entry_id}.{ext}"
        r = requests.get(url, timeout=30)
        if r.status_code == 404:
            continue
        r.raise_for_status()
        dest.write_bytes(r.content)
        return dest
    raise FileNotFoundError(f"{entry_id}: no .pdb or .cif available from RCSB")


def uniprot_protein_family(uniprot_ids: list[str], cfg: dict) -> dict:
    """Batch-fetch UniProt's free-text family membership + broad keywords per
    accession, for stratifying downstream results (e.g. kinases vs. everything
    else) without needing a join at analysis time.

    Family membership comes back under commentType "SIMILARITY" (there is no
    "FAMILY" commentType in the UniProt REST schema — confirmed live against
    rest.uniprot.org for P00915/P00918, both of which return their "Belongs
    to the alpha-carbonic anhydrase family" line under SIMILARITY)."""
    base_url = "https://rest.uniprot.org/uniprotkb/accessions"
    batch_size = cfg["rcsb"]["entry_batch_size"]
    out = {}
    unique_ids = sorted(set(uniprot_ids))
    for i in range(0, len(unique_ids), batch_size):
        batch = unique_ids[i:i + batch_size]
        params = {
            "accessions": ",".join(batch),
            "fields": "accession,protein_families,keyword",
        }
        r = requests.get(base_url, params=params, timeout=30)
        r.raise_for_status()
        for entry in r.json().get("results", []):
            acc = entry.get("primaryAccession")
            family_comments = [
                c.get("texts", [{}])[0].get("value")
                for c in entry.get("comments", [])
                if c.get("commentType") == "SIMILARITY"
            ]
            keywords = [
                kw.get("name")
                for kw in entry.get("keywords", [])
                if kw.get("category") in ("Molecular function", "Biological process")
            ]
            out[acc] = {
                "protein_family": family_comments[0] if family_comments else None,
                "keywords": ";".join(keywords),
            }
        time.sleep(cfg["rcsb"]["request_pause_s"])
    return out


# ---------------------------------------------------------------------------
# DoGSite3 (proteins.plus) — cross-validation, schema not yet confirmed live.
# ---------------------------------------------------------------------------

def run_dogsite3(entry_id: str, cfg: dict) -> dict:
    """Submit `entry_id` to the DoGSite3 REST API and poll for the result.

    NOT YET WORKING: every payload shape tried against
    https://proteins.plus/api/dogsite3_rest returned
    '400 Invalid number of parameters or incorrect parameter name'.
    Get a working example (proteins.plus issues an API key with one) from
    https://proteins.plus/api and fix the `payload=` below before use.
    Raises NotImplementedError until then so callers fail loudly.
    """
    raise NotImplementedError(
        "DoGSite3 REST payload schema unconfirmed — see docstring / config.yaml [dogsite3]"
    )
