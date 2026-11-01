#!/usr/bin/env python
"""
Phase 2, step 2 - purifying selection, contrasting UL23 with UL30 (Aim 2).

H2: UL23 thymidine kinase is dispensable for replication in culture and is the
primary target of aciclovir activation, so it should show weaker purifying
selection than UL30, which encodes the essential DNA polymerase.

Method. Each isolate CDS is compared to the reference CDS codon by codon, using
the protein alignment to set the reading frame and reference coordinates. Each
codon pair contributes synonymous and nonsynonymous differences counted by the
Nei-Gojobori method with the standard Jukes-Cantor correction, against the
synonymous and nonsynonymous site counts of the reference codon.

This is pairwise-to-reference rather than a full phylogenetic estimate. Within a
species that is the appropriate scale: the quantity of interest is polymorphism
within a population, pN/pS, not divergence between species, dN/dS. Values are
therefore reported as pN/pS throughout and must not be read as dN/dS.

Significance is assessed by bootstrapping over isolates, which respects the fact
that codons within one isolate are not independent observations.

Outputs:
  data/processed/selection_by_gene.tsv
  data/processed/selection_by_site.tsv
  reports/phase2/selection_summary.json
"""

import csv
import json
import math
import os
import random
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Data import CodonTable

ROOT = Path(__file__).resolve().parent.parent
GENE_DIR = ROOT / "data" / "processed" / "genes"
REF_DIR = ROOT / "data" / "raw" / "reference"
PROC = ROOT / "data" / "processed"
REP = ROOT / "reports" / "phase2"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
REFSEQ = {"HSV1": "NC_001806", "HSV2": "NC_001798"}
N_BOOTSTRAP = 1000
SEED = 20261102
# pN/pS inherits the sampling error of its denominator. Below this many pooled
# synonymous differences the ratio is dominated by a handful of segregating
# silent sites and is not interpretable, however tight the bootstrap looks.
MIN_SYN_DIFFS = 50

TABLE = CodonTable.unambiguous_dna_by_id[1].forward_table
BASES = "TCAG"


def translate(codon):
    if len(codon) != 3 or any(b not in BASES for b in codon):
        return None
    return TABLE.get(codon, "*")


def site_counts(codon):
    """Nei-Gojobori synonymous (S) and nonsynonymous (N) site counts."""
    aa = translate(codon)
    if aa is None:
        return None
    syn = 0.0
    for i in range(3):
        for b in BASES:
            if b == codon[i]:
                continue
            alt = codon[:i] + b + codon[i + 1:]
            if translate(alt) == aa:
                syn += 1 / 3
    return syn, 3 - syn


def pathway_diffs(c1, c2):
    """Synonymous and nonsynonymous differences between two codons, averaged
    over all shortest mutational pathways."""
    diff = [i for i in range(3) if c1[i] != c2[i]]
    if not diff:
        return 0.0, 0.0
    if translate(c1) is None or translate(c2) is None:
        return None
    if len(diff) == 1:
        return ((1.0, 0.0) if translate(c1) == translate(c2)
                else (0.0, 1.0))

    import itertools
    sd = nd = 0.0
    n_paths = 0
    for order in itertools.permutations(diff):
        cur, s, n, ok = c1, 0.0, 0.0, True
        for i in order:
            nxt = cur[:i] + c2[i] + cur[i + 1:]
            a, b = translate(cur), translate(nxt)
            if a is None or b is None or a == "*" or b == "*":
                ok = False
                break
            if a == b:
                s += 1
            else:
                n += 1
            cur = nxt
        if ok:
            sd += s
            nd += n
            n_paths += 1
    if not n_paths:
        return None
    return sd / n_paths, nd / n_paths


def jukes_cantor(p):
    """JC correction; returns None where the proportion is out of range."""
    if p is None or p <= 0:
        return 0.0
    if p >= 0.75:
        return None
    return -0.75 * math.log(1 - 4 * p / 3)


def read_fasta(path):
    name, buf, out = None, [], []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def reference_cds():
    """Reference nucleotide CDS per (virus, gene), cached on disk."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from p1_02_gate2 import match_gene, parse_cds_fasta, qualifiers

    out = {}
    for virus, acc in REFSEQ.items():
        cache = REF_DIR / f"{virus}_cds_na.fasta"
        if cache.exists() and cache.stat().st_size > 0:
            text = cache.read_text(encoding="utf-8")
        else:
            params = {"db": "nuccore", "id": acc, "rettype": "fasta_cds_na",
                      "retmode": "text", "tool": "hsv-resistance-protocol"}
            if os.environ.get("NCBI_API_KEY"):
                params["api_key"] = os.environ["NCBI_API_KEY"]
            url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(params)
            with urllib.request.urlopen(url, timeout=180) as r:
                text = r.read().decode("utf-8", errors="replace")
            cache.write_text(text, encoding="utf-8")
            time.sleep(0.4)
        for header, seq in parse_cds_fasta(text):
            g = match_gene(qualifiers(header))
            if g and (virus, g) not in out:
                out[(virus, g)] = seq.upper()
    return out


def make_aligner():
    a = Align.PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = "global"
    a.end_insertion_score = 0.0
    a.end_deletion_score = 0.0
    return a


def analyse_gene(virus, gene, ref_nt, aligner):
    """Per-isolate Sd/Nd/S/N totals, plus per-site nonsynonymous counts.

    Differences are counted against the POPULATION CONSENSUS codon at each
    reference position, not against the reference strain itself.

    This matters. pN/pS is a statistic about polymorphism within a population,
    so the comparison point must represent that population. Measuring against
    the reference strain instead counts every position where the reference
    happens to carry a minority allele as a difference in almost every isolate.
    HSV-2 UL30 has three such positions, each present in about 93% of isolates
    and each nonsynonymous; against the HG52 reference they contributed roughly
    510 nonsynonymous differences across 183 isolates and pushed pN/pS above 1
    for an essential polymerase, which is not biologically credible. Against
    the consensus they correctly contribute nothing, because they are not
    polymorphic.

    The reference is still used to define coordinates and the reading frame.
    """
    faa = GENE_DIR / f"{virus}_{gene}.faa"
    fna = GENE_DIR / f"{virus}_{gene}.fna"
    if not (faa.exists() and fna.exists()):
        return None

    prots = dict(read_fasta(faa))
    nts = dict(read_fasta(fna))
    ref_aa = "".join((REF_DIR / f"{virus}_{gene}.faa")
                     .read_text(encoding="utf-8").splitlines()[1:])
    n_ref_codons = len(ref_nt) // 3

    # Pass 1 - map every isolate onto reference coordinates, once.
    mapped = {}
    codon_pool = defaultdict(Counter)
    for acc, aa in prots.items():
        nt = nts.get(acc)
        if not nt:
            continue
        codons = [nt[i:i + 3] for i in range(0, len(nt) - 2, 3)]
        if len(codons) < len(aa):
            continue
        aln = aligner.align(ref_aa, aa)[0]
        m = {}
        for (rs, re_), (qs, _) in zip(*aln.aligned):
            for off in range(re_ - rs):
                ri, qi = rs + off, qs + off
                if ri >= n_ref_codons or qi >= len(codons):
                    continue
                c = codons[qi]
                if translate(c) is None:
                    continue
                m[ri] = c
                codon_pool[ri][c] += 1
        if m:
            mapped[acc] = m

    # Consensus codon per reference position.
    consensus = {ri: pool.most_common(1)[0][0]
                 for ri, pool in codon_pool.items()}
    per_codon_sites = {ri: site_counts(c) for ri, c in consensus.items()}
    S_ref = sum(s[0] for s in per_codon_sites.values() if s)
    N_ref = sum(s[1] for s in per_codon_sites.values() if s)

    # How far the reference strain itself departs from its own population.
    ref_codons = [ref_nt[i:i + 3] for i in range(0, len(ref_nt) - 2, 3)]
    ref_minority = sorted(
        ri + 1 for ri, c in consensus.items()
        if ri < len(ref_codons) and translate(ref_codons[ri]) != translate(c))

    # Pass 2 - differences against the consensus.
    per_isolate = []
    site_nonsyn, site_syn, site_cov = (defaultdict(int), defaultdict(int),
                                       defaultdict(int))
    for acc, m in mapped.items():
        Sd = Nd = S = N = 0.0
        for ri, c2 in m.items():
            sc = per_codon_sites.get(ri)
            if sc is None:
                continue
            c1 = consensus[ri]
            S += sc[0]
            N += sc[1]
            site_cov[ri + 1] += 1
            d = pathway_diffs(c1, c2)
            if d is None:
                continue
            Sd += d[0]
            Nd += d[1]
            if d[1] > 0:
                site_nonsyn[ri + 1] += 1
            if d[0] > 0:
                site_syn[ri + 1] += 1
        if S > 0 and N > 0:
            per_isolate.append({"accession": acc, "Sd": Sd, "Nd": Nd,
                                "S": S, "N": N})

    return {
        "per_isolate": per_isolate,
        "site_nonsyn": dict(site_nonsyn),
        "site_syn": dict(site_syn),
        "site_cov": dict(site_cov),
        "S_ref": S_ref, "N_ref": N_ref,
        "n_codons": n_ref_codons,
        "n_reference_minority_positions": len(ref_minority),
        "reference_minority_positions": ref_minority[:50],
    }


def pn_ps(per_isolate):
    """Pooled pN, pS and their ratio, with Jukes-Cantor correction."""
    Sd = sum(x["Sd"] for x in per_isolate)
    Nd = sum(x["Nd"] for x in per_isolate)
    S = sum(x["S"] for x in per_isolate)
    N = sum(x["N"] for x in per_isolate)
    if S == 0 or N == 0:
        return None
    pS = jukes_cantor(Sd / S)
    pN = jukes_cantor(Nd / N)
    if pS is None or pN is None or pS == 0:
        return {"pN": pN, "pS": pS, "ratio": None,
                "Sd": Sd, "Nd": Nd, "S": S, "N": N}
    return {"pN": pN, "pS": pS, "ratio": pN / pS,
            "Sd": Sd, "Nd": Nd, "S": S, "N": N}


def bootstrap_ratio(per_isolate, n=N_BOOTSTRAP, seed=SEED):
    rng = random.Random(seed)
    vals = []
    k = len(per_isolate)
    if k < 10:
        return None
    for _ in range(n):
        sample = [per_isolate[rng.randrange(k)] for _ in range(k)]
        r = pn_ps(sample)
        if r and r.get("ratio") is not None:
            vals.append(r["ratio"])
    if not vals:
        return None
    vals.sort()
    return {"lo": vals[int(0.025 * len(vals))],
            "hi": vals[int(0.975 * len(vals)) - 1],
            "n_resamples": len(vals)}


def main():
    REP.mkdir(parents=True, exist_ok=True)
    aligner = make_aligner()
    refs_nt = reference_cds()

    gene_rows, site_rows, summary = [], [], {}
    store = {}

    for (virus, gene), ref_nt in sorted(refs_nt.items()):
        res = analyse_gene(virus, gene, ref_nt, aligner)
        if not res or not res["per_isolate"]:
            continue
        est = pn_ps(res["per_isolate"])
        boot = bootstrap_ratio(res["per_isolate"])
        store[(virus, gene)] = res

        # Concentration of the nonsynonymous signal: if a few sites carry most
        # of it, the gene-wide ratio describes those sites, not the gene.
        ns = sorted(res["site_nonsyn"].values(), reverse=True)
        top5_share = (sum(ns[:5]) / sum(ns)) if ns and sum(ns) else 0.0
        reliable = est["Sd"] >= MIN_SYN_DIFFS and top5_share < 0.6

        gene_rows.append({
            "virus": virus, "gene": gene,
            "n_isolates": len(res["per_isolate"]),
            "n_codons": res["n_codons"],
            "syn_sites": round(res["S_ref"], 1),
            "nonsyn_sites": round(res["N_ref"], 1),
            "Sd": round(est["Sd"], 2), "Nd": round(est["Nd"], 2),
            "pS": round(est["pS"], 6) if est["pS"] is not None else "",
            "pN": round(est["pN"], 6) if est["pN"] is not None else "",
            "pN_pS": round(est["ratio"], 4) if est.get("ratio") else "",
            "ci_low": round(boot["lo"], 4) if boot else "",
            "ci_high": round(boot["hi"], 4) if boot else "",
            "top5_site_share_of_nonsyn": round(top5_share, 3),
            "reliable": int(reliable),
        })
        summary[f"{virus}:{gene}"] = {
            "pN_pS": est.get("ratio"),
            "ci": [boot["lo"], boot["hi"]] if boot else None,
            "n_isolates": len(res["per_isolate"]),
            "Sd": round(est["Sd"], 1),
            "top5_site_share_of_nonsyn": round(top5_share, 3),
            "reliable": reliable,
            "reference_minority_positions":
                res["n_reference_minority_positions"],
        }
        print(f"{virus}:{gene:<6} n={len(res['per_isolate']):<4} "
              f"pN/pS={est.get('ratio'):.4f} "
              f"[{boot['lo']:.4f}, {boot['hi']:.4f}]" if boot and est.get("ratio")
              else f"{virus}:{gene:<6} pN/pS unavailable")

        for pos in sorted(res["site_cov"]):
            cov = res["site_cov"][pos]
            if cov == 0:
                continue
            site_rows.append({
                "virus": virus, "gene": gene, "position": pos,
                "coverage": cov,
                "n_nonsyn": res["site_nonsyn"].get(pos, 0),
                "n_syn": res["site_syn"].get(pos, 0),
                "nonsyn_frac": round(res["site_nonsyn"].get(pos, 0) / cov, 5),
            })

    # The Aim 2 contrast: dispensable kinase against essential polymerase.
    contrasts = {}
    for virus in REFSEQ:
        a = summary.get(f"{virus}:UL23")
        b = summary.get(f"{virus}:UL30")
        if not (a and b and a["pN_pS"] and b["pN_pS"]):
            continue
        if not (a["reliable"] and b["reliable"]):
            contrasts[virus] = {"skipped": "one or both genes unreliable",
                                "UL23_reliable": a["reliable"],
                                "UL30_reliable": b["reliable"]}
            continue
        res_a = store[(virus, "UL23")]["per_isolate"]
        res_b = store[(virus, "UL30")]["per_isolate"]
        rng = random.Random(SEED)
        diffs = []
        for _ in range(N_BOOTSTRAP):
            sa = [res_a[rng.randrange(len(res_a))] for _ in range(len(res_a))]
            sb = [res_b[rng.randrange(len(res_b))] for _ in range(len(res_b))]
            ra, rb = pn_ps(sa), pn_ps(sb)
            if ra and rb and ra.get("ratio") and rb.get("ratio"):
                diffs.append(ra["ratio"] - rb["ratio"])
        diffs.sort()
        contrasts[virus] = {
            "UL23_pN_pS": a["pN_pS"], "UL30_pN_pS": b["pN_pS"],
            "difference": a["pN_pS"] - b["pN_pS"],
            "ci_low": diffs[int(0.025 * len(diffs))] if diffs else None,
            "ci_high": diffs[int(0.975 * len(diffs)) - 1] if diffs else None,
            "fraction_of_resamples_UL23_higher":
                round(sum(1 for d in diffs if d > 0) / len(diffs), 4)
                if diffs else None,
        }

    for name, data in (("selection_by_gene.tsv", gene_rows),
                       ("selection_by_site.tsv", site_rows)):
        if not data:
            continue
        with open(PROC / name, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0].keys()),
                               delimiter="\t")
            w.writeheader()
            w.writerows(data)

    (REP / "selection_summary.json").write_text(
        json.dumps({"per_gene": summary, "aim2_contrast": contrasts,
                    "n_bootstrap": N_BOOTSTRAP, "seed": SEED},
                   indent=2, sort_keys=True), encoding="utf-8")

    print("\nAim 2 contrast, UL23 against UL30:")
    for virus, c in sorted(contrasts.items()):
        if "skipped" in c:
            print(f"  {virus}: not computed - {c['skipped']}")
            continue
        print(f"  {virus}: UL23 {c['UL23_pN_pS']:.4f} vs UL30 "
              f"{c['UL30_pN_pS']:.4f}  difference {c['difference']:+.4f} "
              f"[{c['ci_low']:+.4f}, {c['ci_high']:+.4f}]  "
              f"support {c['fraction_of_resamples_UL23_higher']:.1%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
