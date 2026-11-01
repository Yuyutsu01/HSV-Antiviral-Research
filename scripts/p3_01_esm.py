#!/usr/bin/env python
"""
Phase 3, step 1 - ESM-2 masked-marginal scores for every position in the four
target proteins.

Used zero-shot: no fine-tuning, no training on resistance labels. The score for
substituting reference residue r with alternate a at position i is

    log P(a | context, position i masked) - log P(r | context, position i masked)

which is the standard masked-marginal formulation. A more negative score means
the model considers the substitution less plausible given the rest of the
protein.

This feature is the main circularity risk in the whole design: a protein
language model learns the same evolutionary constraint that the conservation
features measure directly, so any apparent contribution has to be checked
against conservation before it is believed. The audit in p3_05 does that.

Scores every position of every target protein once and caches the result, so
the model runs even though only ~650 variants are ultimately scored.

Outputs:
  data/processed/esm2_scores.tsv
"""

import csv
import sys
import time
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parent.parent
REF_DIR = ROOT / "data" / "raw" / "reference"
OUT = ROOT / "data" / "processed"

MODEL_NAME = "esm2_t33_650M_UR50D"
AA = "ACDEFGHIKLMNPQRSTVWY"
# Long proteins are scored in overlapping windows; ESM-2 handles 1022 tokens.
WINDOW = 1000
STRIDE = 500
# A batch of 32 windows of 1000 tokens overflows a 6 GB card and spills into
# host memory, which costs an order of magnitude in speed. 8 stays resident.
BATCH = 8
# Scoring is one forward pass per position, so each protein is written as it
# finishes and a restart resumes from there rather than repeating the work.
CHECKPOINT_EVERY_PROTEIN = True


def read_fasta_one(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    return "".join(l.strip() for l in lines[1:] if l.strip())


def windows(n):
    """Start offsets covering 0..n, with each position scored in the window
    where it sits furthest from an edge."""
    if n <= WINDOW:
        return [(0, n)]
    out = []
    start = 0
    while start < n:
        out.append((start, min(start + WINDOW, n)))
        if start + WINDOW >= n:
            break
        start += STRIDE
    return out


FIELDS = ["virus", "gene", "position", "ref_aa", "alt_aa",
          "esm2_masked_marginal", "esm2_ref_logprob"]


def load_done(path):
    """Proteins already present in a partial output, so a restart resumes."""
    if not path.exists():
        return set(), []
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {(r["virus"], r["gene"]) for r in rows}, rows


def write_rows(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"loading {MODEL_NAME} on {device} (first run downloads ~2.5 GB)",
          flush=True)
    import esm
    model, alphabet = esm.pretrained.__dict__[MODEL_NAME]()
    model = model.to(device).eval()
    if device == "cuda":
        model = model.half()          # halves activation memory, keeps it resident
    bc = alphabet.get_batch_converter()
    aa_idx = {a: alphabet.get_idx(a) for a in AA}

    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / "esm2_scores.tsv"
    done, rows = load_done(out_path)
    if done:
        print(f"resuming, already scored: {sorted(done)}", flush=True)

    targets = sorted(p for p in REF_DIR.glob("*_*.faa")
                     if not p.stem.endswith("_cds_aa"))
    for path in targets:
        virus, gene = path.stem.split("_", 1)
        if (virus, gene) in done:
            continue
        seq = read_fasta_one(path)
        n = len(seq)
        t0 = time.time()
        # best[i] = (distance_from_window_edge, logprob_vector)
        best = {}
        for (s, e) in windows(n):
            sub = seq[s:e]
            _, _, toks = bc([("x", sub)])
            toks = toks.to(device)
            # Mask one position at a time, batched for speed.
            for batch_start in range(0, len(sub), BATCH):
                idxs = list(range(batch_start,
                                  min(batch_start + BATCH, len(sub))))
                rep = toks.repeat(len(idxs), 1)
                for r, i in enumerate(idxs):
                    rep[r, i + 1] = alphabet.mask_idx   # +1 for BOS
                with torch.no_grad():
                    logits = model(rep)["logits"]
                # log_softmax in float32 regardless of the model's precision.
                lp = torch.log_softmax(logits.float(), dim=-1)
                for r, i in enumerate(idxs):
                    pos = s + i
                    edge = min(i, len(sub) - 1 - i)
                    if pos in best and best[pos][0] >= edge:
                        continue
                    best[pos] = (edge, lp[r, i + 1].detach().cpu())
                if batch_start % (BATCH * 20) == 0:
                    frac = (batch_start + len(idxs)) / len(sub)
                    print(f"  {virus}:{gene} window {s}-{e} "
                          f"{frac:5.1%} ({time.time() - t0:.0f}s)", flush=True)
            print(f"  {virus}:{gene} window {s}-{e} done", flush=True)

        for pos in range(n):
            ref_aa = seq[pos]
            if ref_aa not in aa_idx:
                continue
            lp = best[pos][1]
            ref_lp = lp[aa_idx[ref_aa]].item()
            for alt in AA:
                if alt == ref_aa:
                    continue
                rows.append({
                    "virus": virus, "gene": gene, "position": pos + 1,
                    "ref_aa": ref_aa, "alt_aa": alt,
                    "esm2_masked_marginal": round(
                        lp[aa_idx[alt]].item() - ref_lp, 5),
                    "esm2_ref_logprob": round(ref_lp, 5),
                })
        print(f"{virus}:{gene} scored, {n} positions "
              f"in {time.time() - t0:.0f}s", flush=True)
        if CHECKPOINT_EVERY_PROTEIN:
            write_rows(out_path, rows)

    write_rows(out_path, rows)
    print(f"\n{len(rows)} substitution scores written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
