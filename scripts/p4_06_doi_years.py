"""Phase 4, step 6 — resolve HerpesDRG reference DOIs to publication years.

Needed for the temporal hold-out, which the protocol calls the most persuasive
validation in the plan but which could not be run: publication year is
recoverable from the DOI string itself for only 35 of 89 distinct studies, and
the recoverable ones cluster by publisher, so a split on partial years would
confound era with journal.

This queries Crossref for the remainder. It writes a cache with its own
retrieval date, so the lookup runs once and the temporal analysis is afterwards
reproducible offline. Re-run with --force to refresh.

Only DOIs are sent. No mailto parameter and no address in the User-Agent: that
would put a personal identifier into third-party request logs for no analytical
gain. The cost is the anonymous rate-limit pool, which for 89 requests is
irrelevant.

The year taken is the earliest of `published-print` and `published-online`,
because a 2015 online-first paper that appeared in a 2016 issue was knowable to
the field in 2015, and the hold-out asks what was knowable by a date.
"""

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB = os.path.join(ROOT, "data", "raw", "herpesdrg-db.tsv")
CACHE = os.path.join(ROOT, "data", "raw", "doi_years.json")

API = "https://api.crossref.org/works/"
UA = "hsv-resistance-frequencies/1.0 (academic research; no contact in UA by design)"
DELAY = 0.2


def normalise(doi):
    d = doi.strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi.org/",
                   "https://dx.doi.org/", "doi:"):
        if d.lower().startswith(prefix):
            d = d[len(prefix):]
    return d.strip().rstrip(".")


def distinct_dois():
    import csv
    with open(DB, encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["status"] == "A" and r["virus"] in ("HSV1", "HSV2")]
    return sorted({normalise(r["ref_doi"]) for r in rows if r["ref_doi"].strip()})


def year_from_message(msg):
    """Earliest of published-print and published-online."""
    years = []
    for field in ("published-online", "published-print", "published",
                  "issued", "created"):
        part = msg.get(field) or {}
        dp = part.get("date-parts") or []
        if dp and dp[0] and isinstance(dp[0][0], int):
            years.append((dp[0][0], field))
    if not years:
        return None, None
    years.sort()
    return years[0]


def fetch(doi, tries=3):
    req = urllib.request.Request(API + urllib.parse.quote(doi, safe="/"),
                                 headers={"User-Agent": UA,
                                          "Accept": "application/json"})
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode("utf-8"))["message"]
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code in (429, 500, 502, 503, 504) and attempt < tries - 1:
                time.sleep(2 ** attempt)
                continue
            raise
        except (urllib.error.URLError, TimeoutError):
            if attempt < tries - 1:
                time.sleep(2 ** attempt)
                continue
            raise
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    dois = distinct_dois()
    print(f"{len(dois)} distinct DOIs across active HSV rows", flush=True)

    cache = {"retrieved": None, "source": "api.crossref.org", "years": {}}
    if os.path.exists(CACHE) and not args.force:
        with open(CACHE, encoding="utf-8") as fh:
            cache = json.load(fh)
    known = cache["years"]

    todo = [d for d in dois if d not in known]
    print(f"{len(known)} cached, {len(todo)} to fetch", flush=True)

    failures = []
    for i, doi in enumerate(todo, 1):
        try:
            msg = fetch(doi)
        except Exception as exc:                                # noqa: BLE001
            failures.append((doi, repr(exc)))
            print(f"  [{i}/{len(todo)}] {doi}  ERROR {exc}", flush=True)
            time.sleep(DELAY)
            continue
        if msg is None:
            known[doi] = {"year": None, "note": "not found in Crossref"}
            print(f"  [{i}/{len(todo)}] {doi}  not found", flush=True)
        else:
            year, field = year_from_message(msg)
            known[doi] = {
                "year": year,
                "year_field": field,
                "container": (msg.get("container-title") or [""])[0],
                "publisher": msg.get("publisher", ""),
                "title": (msg.get("title") or [""])[0][:160],
            }
            print(f"  [{i}/{len(todo)}] {doi}  {year}  "
                  f"{known[doi]['container'][:48]}", flush=True)
        time.sleep(DELAY)

    cache["retrieved"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    cache["years"] = known
    with open(CACHE, "w", encoding="utf-8") as fh:
        json.dump(cache, fh, indent=2, sort_keys=True)

    resolved = {d: v["year"] for d, v in known.items()
                if d in dois and v.get("year")}
    print(f"\nresolved {len(resolved)}/{len(dois)} DOIs to a year", flush=True)
    if failures:
        print(f"{len(failures)} request failure(s):", flush=True)
        for d, e in failures:
            print(f"   {d}  {e}", flush=True)
    missing = [d for d in dois if d not in resolved]
    if missing:
        print(f"{len(missing)} unresolved:", flush=True)
        for d in missing:
            print(f"   {d}  {known.get(d, {}).get('note', '')}", flush=True)

    print("\nyear distribution:", flush=True)
    for y, n in sorted(Counter(resolved.values()).items()):
        print(f"   {y}  {'#' * n} {n}", flush=True)
    print(f"\npublishers represented: "
          f"{len({v.get('publisher') for v in known.values() if v.get('year')})}",
          flush=True)
    print(f"wrote {CACHE}", flush=True)
    return 0 if len(resolved) >= 0.9 * len(dois) else 1


if __name__ == "__main__":
    sys.exit(main())
