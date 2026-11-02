"""Phase 4, step 3 — the release artefacts.

The protocol's deposition list: a frozen accession list with its query date, the
master variant table, and the reclassification table deposited as a reusable
resource rather than buried in a supplement.

Three things this script is careful about.

**The accession list records what was excluded, not only what was kept.** A
frozen list of the 451 analysed isolates is not reproducible on its own: a
reader cannot tell whether an accession is absent because the query missed it
or because a filter dropped it. Every accession returned by the original query
appears here with its fate.

**The master variant table is the union of two populations, not their
intersection.** Catalogued mutations and observed population variants overlap
only partly -- that disjunction is the project's central finding -- so a table
restricted to either one would misrepresent it. `in_catalogue` and `observed`
are separate columns and either can be zero.

**Features are carried for every row that has them**, independently of whether
the variant has a phenotype. Only 235 variants have a measured fold change,
but conservation, structure and ESM-2 features exist for far more, and a
reader reproducing Aim 4 needs to see the rows that were available but
unlabelled.
"""

import csv
import hashlib
import json
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")
RAW = os.path.join(ROOT, "data", "raw")
INTERIM = os.path.join(ROOT, "data", "interim")
REL = os.path.join(ROOT, "data", "release")
REP = os.path.join(ROOT, "reports", "phase4")


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------- accessions
def build_accession_list():
    queries, manifest = {}, {}
    for virus in ("HSV1", "HSV2"):
        with open(os.path.join(RAW, "accessions", f"{virus}.json"),
                  encoding="utf-8") as fh:
            d = json.load(fh)
        queries[virus] = {"taxid": d["taxid"], "query_term": d["query_term"],
                          "query_date": d["query_date"], "n_returned": d["n"]}
        for acc in d["accessions"]:
            manifest[acc] = virus

    meta = {r["accession"]: r for r in read_tsv(
        os.path.join(INTERIM, "isolate_metadata.tsv"))}

    analysed = defaultdict(set)
    for r in read_tsv(os.path.join(PROC, "isolate_coverage.tsv")):
        analysed[r["accession"]].add(r["gene"])

    collapsed = {}
    cpath = os.path.join(INTERIM, "patient_collapse.tsv")
    if os.path.exists(cpath):
        for r in read_tsv(cpath):
            acc = r.get("accession") or r.get("dropped_accession") or ""
            if acc:
                collapsed[acc] = r

    # Per-gene CDS QC, so a failure can name its reason instead of falling into
    # a catch-all bucket. An accession that never reached CDS extraction is a
    # different failure from one whose extracted genes failed QC.
    qc_flags, qc_any_pass = defaultdict(set), defaultdict(bool)
    for r in read_tsv(os.path.join(INTERIM, "cds_qc.tsv")):
        if r["passed"] == "1":
            qc_any_pass[r["accession"]] = True
        else:
            for fl in r["flags"].split(";"):
                if fl and fl != "pass":
                    qc_flags[r["accession"]].add(fl)

    rows = []
    for acc, virus in sorted(manifest.items()):
        m = meta.get(acc, {})
        genes = analysed.get(acc, set())
        if genes:
            fate, detail = "analysed", ""
        elif acc in collapsed:
            fate, detail = "excluded_patient_duplicate", collapsed[acc].get(
                "kept_representative", "")
        elif str(m.get("is_lab_strain", "")) in ("1", "True", "true"):
            fate, detail = "excluded_laboratory_strain", m.get("strain", "")
        elif m and m.get("host") and "homo sapiens" not in m["host"].lower():
            fate, detail = "excluded_non_human_host", m.get("host", "")
        elif acc in qc_flags:
            fate = "excluded_cds_qc"
            detail = ";".join(sorted(qc_flags[acc]))
        elif acc not in qc_any_pass:
            fate, detail = "excluded_no_cds_extracted", ""
        elif m:
            fate, detail = "excluded_downstream_filter", ""
        else:
            fate, detail = "excluded_no_metadata", ""
        rows.append({
            "accession": acc, "virus": virus, "fate": fate,
            "fate_detail": detail,
            "n_genes_analysed": len(genes),
            "genes_analysed": ";".join(sorted(genes)),
            "study_key": m.get("study_key", ""),
            "pubmed": m.get("pubmed", ""),
            "collection_year": m.get("collection_year", ""),
            "country": m.get("country", ""),
            "host": m.get("host", ""),
            "is_lab_strain": m.get("is_lab_strain", ""),
            "strain": m.get("strain", ""),
            "query_term": queries[virus]["query_term"],
            "query_date": queries[virus]["query_date"],
        })

    path = write_tsv(
        os.path.join(REL, "accession_list.tsv"),
        ["accession", "virus", "fate", "fate_detail", "n_genes_analysed",
         "genes_analysed",
         "study_key", "pubmed", "collection_year", "country", "host",
         "is_lab_strain", "strain", "query_term", "query_date"], rows)

    fates = defaultdict(int)
    for r in rows:
        fates[r["fate"]] += 1
    return path, queries, dict(fates), len(rows)


# ----------------------------------------------------------- master variants
def build_master_table():
    key = lambda r: (r["virus"], r["gene"], r["position"],
                     r["ref_aa"], r["alt_aa"])

    cat = {key(r): r for r in read_tsv(
        os.path.join(PROC, "reclassification_table.tsv"))}
    freq = {key(r): r for r in read_tsv(
        os.path.join(PROC, "allele_frequencies.tsv"))}
    esm = {key(r): r for r in read_tsv(os.path.join(PROC, "esm2_scores.tsv"))}
    struct = {(r["virus"], r["gene"], r["position"]): r for r in read_tsv(
        os.path.join(PROC, "structure_features.tsv"))}
    site = {(r["virus"], r["gene"], r["position"]): r for r in read_tsv(
        os.path.join(PROC, "selection_by_site.tsv"))}
    model = {(r["virus"], r["gene"], r["position"], r["ref_aa"], r["alt_aa"]): r
             for r in read_tsv(os.path.join(PROC, "model_dataset.tsv"))}

    rows = []
    for k in sorted(set(cat) | set(freq)):
        virus, gene, pos, ref, alt = k
        c, f = cat.get(k, {}), freq.get(k, {})
        s = struct.get((virus, gene, pos), {})
        t = site.get((virus, gene, pos), {})
        e = esm.get(k, {})
        m = model.get(k, {})
        rows.append({
            "virus": virus, "gene": gene, "position": pos,
            "ref_aa": ref, "alt_aa": alt,
            "aa_change": c.get("aa_change") or f"{ref}{pos}{alt}",
            # --- catalogue side
            "in_catalogue": int(bool(c)),
            "evidence_tier": c.get("evidence_tier", ""),
            "n_refs": c.get("n_refs", ""),
            "drugs_tested": c.get("drugs_tested", ""),
            "max_fold_change": c.get("max_fold_change", ""),
            "resistant_at_2x": m.get("resistant", ""),
            # --- population side
            "observed": int(bool(f)) if f else int(c.get("observed") or 0),
            "count": f.get("count", c.get("count", "")),
            "coverage": f.get("coverage", c.get("coverage", "")),
            "frequency": f.get("frequency", c.get("frequency", "")),
            "ci_low": f.get("ci_low", c.get("ci_low", "")),
            "ci_high": f.get("ci_high", c.get("ci_high", "")),
            # --- ascertainment, from the reclassification analysis
            "n_studies_carrying": c.get("n_studies_carrying", ""),
            "loso_min_freq": c.get("loso_min_freq", ""),
            "max_study_share": c.get("max_study_share", ""),
            "study_dominated": c.get("study_dominated", ""),
            "frequency_band": c.get("frequency_band", ""),
            "flag": c.get("flag", ""),
            # --- features
            "cons_site_variability": t.get("nonsyn_frac", ""),
            "cons_site_coverage": t.get("coverage", ""),
            "struct_rsa": s.get("rsa", ""),
            "struct_dist_to_ligand": s.get("dist_to_ligand", ""),
            "struct_site_resolved": s.get("site_resolved", ""),
            "struct_pdb_id": s.get("pdb_id", ""),
            "struct_inherited": s.get("inherited", ""),
            "esm2_masked_marginal": e.get("esm2_masked_marginal", ""),
            "esm2_ref_logprob": e.get("esm2_ref_logprob", ""),
            "in_model_dataset": int(bool(m)),
        })

    fields = list(rows[0].keys())
    path = write_tsv(os.path.join(REL, "master_variant_table.tsv"), fields, rows)
    stats = {
        "n_rows": len(rows),
        "n_in_catalogue": sum(r["in_catalogue"] for r in rows),
        "n_observed": sum(r["observed"] for r in rows),
        "n_both": sum(1 for r in rows if r["in_catalogue"] and r["observed"]),
        "n_catalogue_only": sum(1 for r in rows
                                if r["in_catalogue"] and not r["observed"]),
        "n_observed_only": sum(1 for r in rows
                               if r["observed"] and not r["in_catalogue"]),
        "n_in_model_dataset": sum(r["in_model_dataset"] for r in rows),
        "n_with_structure": sum(1 for r in rows if r["struct_rsa"] != ""),
        "n_with_esm": sum(1 for r in rows if r["esm2_masked_marginal"] != ""),
    }
    return path, stats, fields


def main():
    os.makedirs(REL, exist_ok=True)
    os.makedirs(REP, exist_ok=True)

    acc_path, queries, fates, n_acc = build_accession_list()
    print(f"accession list: {n_acc} accessions", flush=True)
    for k, v in sorted(fates.items(), key=lambda kv: -kv[1]):
        print(f"   {k:34s} {v}", flush=True)

    master_path, stats, fields = build_master_table()
    print(f"\nmaster variant table: {stats['n_rows']} rows, "
          f"{len(fields)} columns", flush=True)
    for k, v in stats.items():
        print(f"   {k:24s} {v}", flush=True)

    # The reclassification table is deposited as its own resource.
    src = os.path.join(PROC, "reclassification_table.tsv")
    dep = os.path.join(REL, "reclassification_table.tsv")
    reclass = read_tsv(src)
    write_tsv(dep, list(reclass[0].keys()), reclass)
    print(f"\ndeposited reclassification table: {len(reclass)} rows", flush=True)

    # ---- gates
    assert stats["n_in_catalogue"] == 644, stats["n_in_catalogue"]
    assert fates.get("analysed", 0) > 0
    assert stats["n_observed"] >= stats["n_both"]

    out = {
        "queries": queries,
        "accession_fates": fates,
        "master_table": stats,
        "files": {},
    }
    for p in (acc_path, master_path, dep):
        rel = os.path.relpath(p, ROOT).replace("\\", "/")
        out["files"][rel] = {"sha256": sha256(p),
                             "bytes": os.path.getsize(p)}
        print(f"   {rel}  sha256={out['files'][rel]['sha256'][:16]}...",
              flush=True)

    path = os.path.join(REP, "release_summary.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
