#!/usr/bin/env python
"""
Phase 1, step 1 - freeze the accession list and fetch sequence data from NCBI.

Two passes per virus:
  A. coding sequences  (efetch rettype=fasta_cds_na) -> data/raw/cds/
  B. record metadata   (efetch rettype=gb, truncated) -> data/raw/meta/

The accession list is written once with its query date and reused on every later
run, so the dataset is frozen at the moment of first download. Delete
data/raw/accessions/*.json to re-query, which starts a new dataset.

Resumable: batches already on disk are skipped. Safe to re-run after an
interruption.

Optional environment variables:
  NCBI_API_KEY   raises the rate limit from 3/s to 10/s
  NCBI_EMAIL     identifies the caller to NCBI, as their usage policy requests
"""

import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ACC_DIR = ROOT / "data" / "raw" / "accessions"
CDS_DIR = ROOT / "data" / "raw" / "cds"
META_DIR = ROOT / "data" / "raw" / "meta"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Complete-genome length window. HSV-1 is ~152 kb and HSV-2 ~155 kb; the window
# is deliberately loose here and tightened per virus during QC in step 2.
VIRUSES = {
    "HSV1": {"taxid": 10298, "refseq": "NC_001806"},
    "HSV2": {"taxid": 10310, "refseq": "NC_001798"},
}
LEN_MIN, LEN_MAX = 140_000, 200_000

API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
EMAIL = os.environ.get("NCBI_EMAIL", "").strip()
DELAY = 0.12 if API_KEY else 0.36          # stay under NCBI's published limits
BATCH_CDS = 10
BATCH_META = 50


def _params(extra):
    p = {"tool": "hsv-resistance-protocol"}
    if API_KEY:
        p["api_key"] = API_KEY
    if EMAIL:
        p["email"] = EMAIL
    p.update(extra)
    return p


def eutil(endpoint, params, retries=5):
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(_params(params))
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"{endpoint} failed after {retries} tries: {last}")


def freeze_accessions(virus):
    """Query once, then reuse. The saved file defines the dataset."""
    ACC_DIR.mkdir(parents=True, exist_ok=True)
    path = ACC_DIR / f"{virus}.json"
    if path.exists():
        d = json.loads(path.read_text(encoding="utf-8"))
        print(f"[{virus}] accession list frozen {d['query_date']}: "
              f"{len(d['accessions'])} records")
        return d["accessions"]

    taxid = VIRUSES[virus]["taxid"]
    term = f"txid{taxid}[Organism:exp] AND {LEN_MIN}:{LEN_MAX}[SLEN]"
    txt = eutil("esearch", {
        "db": "nuccore", "term": term, "retmax": 10000, "retmode": "json",
    })
    ids = json.loads(txt)["esearchresult"]["idlist"]
    time.sleep(DELAY)

    # Resolve UIDs to versioned accessions so the list is human-readable and
    # stable across NCBI's internal id changes.
    accs = []
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        summ = json.loads(eutil("esummary", {
            "db": "nuccore", "id": ",".join(chunk), "retmode": "json",
        }))["result"]
        accs += [summ[u]["accessionversion"] for u in chunk if u in summ]
        time.sleep(DELAY)

    accs = sorted(set(accs))
    path.write_text(json.dumps({
        "virus": virus,
        "taxid": taxid,
        "query_term": term,
        "query_date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n": len(accs),
        "accessions": accs,
    }, indent=2), encoding="utf-8")
    print(f"[{virus}] accession list frozen now: {len(accs)} records")
    return accs


def fetch_batches(virus, accs, outdir, rettype, batch, label, extra=None):
    outdir = outdir / virus
    outdir.mkdir(parents=True, exist_ok=True)
    todo = []
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        f = outdir / f"{i // batch:05d}.txt"
        if not (f.exists() and f.stat().st_size > 0):
            todo.append((f, chunk))
    if not todo:
        print(f"[{virus}] {label}: already complete ({len(accs)} records)")
        return
    print(f"[{virus}] {label}: {len(todo)} batches to fetch")
    for n, (f, chunk) in enumerate(todo, 1):
        params = {"db": "nuccore", "id": ",".join(chunk),
                  "rettype": rettype, "retmode": "text"}
        if extra:
            params.update(extra)
        f.write_text(eutil("efetch", params), encoding="utf-8")
        if n % 10 == 0 or n == len(todo):
            print(f"    {label} {n}/{len(todo)}", flush=True)
        time.sleep(DELAY)


def main():
    for virus in VIRUSES:
        accs = freeze_accessions(virus)
        # Pass A: all coding sequences. Target genes are filtered out in step 3.
        fetch_batches(virus, accs, CDS_DIR, "fasta_cds_na", BATCH_CDS, "CDS")
        # Pass B: record headers only. seq_stop truncates the sequence but keeps
        # the source feature qualifiers (strain, host, country, collection_date).
        fetch_batches(virus, accs, META_DIR, "gb", BATCH_META, "metadata",
                      extra={"seq_start": 1, "seq_stop": 1})
    print("\nfetch complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
