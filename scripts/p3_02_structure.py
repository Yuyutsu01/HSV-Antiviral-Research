#!/usr/bin/env python
"""
Phase 3, step 2 - structural features in reference coordinates.

Two features per position:
  rsa                relative solvent accessibility (Shrake-Rupley, normalised
                     by Tien et al. theoretical maximum ASA per residue)
  dist_to_ligand     minimum heavy-atom distance to any bound ligand

The active site is taken from the ligands actually present in the structure -
the thymidine in a thymidine kinase, the nucleotide and DNA in a polymerase -
rather than from a hand-written list of catalytic residues. That keeps the
feature objective and avoids encoding the published claim about where
resistance mutations sit, which is exactly the circularity this protocol is
trying to measure.

PDB residue numbering routinely differs from the reference protein's, so each
chain is aligned to the reference sequence and features are reported in
reference coordinates. Chains are auto-selected by alignment identity and
coverage across a candidate list.

HSV-2 has no structures for these genes. HSV-2 positions inherit the HSV-1
structure through a pairwise alignment, and the inherited rows are flagged.

Outputs:
  data/processed/structure_features.tsv
  reports/phase3/structure_summary.json
"""

import csv
import json
import sys
import urllib.request
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio.PDB.SASA import ShrakeRupley

ROOT = Path(__file__).resolve().parent.parent
REF_DIR = ROOT / "data" / "raw" / "reference"
PDB_DIR = ROOT / "data" / "raw" / "pdb"
OUT = ROOT / "data" / "processed"
REP = ROOT / "reports" / "phase3"

# Candidates per gene, tried in order; the best-covering chain wins.
CANDIDATES = {
    "UL23": ["1kim", "1vtk", "2ki5"],
    "UL30": ["8exx", "9ja3", "9enq", "8v1r"],
    "UL5": ["9vlq", "9ycp", "9nnp"],
    "UL52": ["9vlq", "9ycp", "9nnp"],
}

# Tien et al. 2013, theoretical maximum ASA.
MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}

# Not active-site ligands: solvent, cryoprotectants and common buffer ions.
IGNORE_HET = {
    "HOH", "DOD", "GOL", "EDO", "PEG", "PG4", "SO4", "PO4", "CL", "NA", "K",
    "MG", "MN", "ZN", "CA", "ACT", "FMT", "TRS", "MES", "EPE", "DMS", "IPA",
}


def fetch_cif(pdb_id):
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    path = PDB_DIR / f"{pdb_id}.cif"
    if path.exists() and path.stat().st_size > 0:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        with urllib.request.urlopen(url, timeout=180) as r:
            path.write_bytes(r.read())
        return path
    except Exception as e:                                # noqa: BLE001
        print(f"    could not fetch {pdb_id}: {e}")
        return None


def read_fasta_one(path):
    return "".join(l.strip() for l in
                   path.read_text(encoding="utf-8").splitlines()[1:] if l.strip())


def make_aligner():
    a = Align.PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = "global"
    a.end_insertion_score = 0.0
    a.end_deletion_score = 0.0
    return a


def chain_sequence(chain):
    """Ordered (residue, one-letter) for standard amino acids in a chain."""
    out = []
    for res in chain:
        if res.id[0] != " ":
            continue
        try:
            aa = protein_letters_3to1[res.get_resname()]
        except KeyError:
            continue
        out.append((res, aa))
    return out


def ligand_atoms(structure):
    out = []
    for res in structure.get_residues():
        het = res.id[0]
        if het == " " or res.get_resname().strip().upper() in IGNORE_HET:
            continue
        out += [a for a in res if a.element != "H"]
    # Nucleic acid chains count as substrate for a polymerase.
    for res in structure.get_residues():
        if res.get_resname().strip() in {"DA", "DT", "DG", "DC", "DU",
                                         "A", "U", "G", "C"}:
            out += [a for a in res if a.element != "H"]
    return out


def analyse(pdb_id, ref_seq, aligner):
    path = fetch_cif(pdb_id)
    if path is None:
        return None
    structure = MMCIFParser(QUIET=True).get_structure(pdb_id, str(path))
    model = next(iter(structure))

    best = None
    for chain in model:
        pairs = chain_sequence(chain)
        if len(pairs) < 50:
            continue
        seq = "".join(a for _, a in pairs)
        aln = aligner.align(ref_seq, seq)[0]
        matched = sum(1 for (rs, re_), (qs, _) in zip(*aln.aligned)
                      for off in range(re_ - rs)
                      if ref_seq[rs + off] == seq[qs + off])
        coverage = matched / len(ref_seq)
        if best is None or coverage > best["coverage"]:
            best = {"chain": chain, "pairs": pairs, "aln": aln,
                    "coverage": coverage, "seq": seq}
    if best is None or best["coverage"] < 0.30:
        return None

    ShrakeRupley().compute(model, level="R")
    lig = ligand_atoms(structure)

    feats = {}
    seq = best["seq"]
    for (rs, re_), (qs, _) in zip(*best["aln"].aligned):
        for off in range(re_ - rs):
            ri, qi = rs + off, qs + off
            if seq[qi] != ref_seq[ri]:
                continue                     # only transfer identical residues
            res = best["pairs"][qi][0]
            aa = seq[qi]
            sasa = getattr(res, "sasa", None)
            rsa = (sasa / MAX_ASA[aa]) if (sasa is not None and aa in MAX_ASA) else None
            d = None
            if lig:
                heavy = [a for a in res if a.element != "H"]
                if heavy:
                    d = min(a - b for a in heavy for b in lig)
            feats[ri + 1] = {"rsa": rsa, "dist_to_ligand": d}
    # In a multi-subunit structure the nearest ligand may belong to another
    # subunit's site, which makes the distance a real number but not an
    # active-site feature for this chain. If no residue of this chain comes
    # within contact range of any ligand, the gene has no resolved site of its
    # own and the distance is flagged as such rather than used as one.
    contacts = sum(1 for f in feats.values()
                   if f["dist_to_ligand"] is not None and f["dist_to_ligand"] < 5.0)
    return {"pdb_id": pdb_id, "chain": best["chain"].id,
            "coverage": round(best["coverage"], 3),
            "n_ligand_atoms": len(lig), "n_contacts": contacts,
            "site_resolved": int(contacts >= 5), "features": feats}


def main():
    REP.mkdir(parents=True, exist_ok=True)
    aligner = make_aligner()
    rows, summary = [], {}

    hsv1_refs = {}
    for gene, cands in CANDIDATES.items():
        ref_path = REF_DIR / f"HSV1_{gene}.faa"
        if not ref_path.exists():
            continue
        ref_seq = read_fasta_one(ref_path)
        hsv1_refs[gene] = ref_seq
        chosen = None
        for pdb_id in cands:
            print(f"HSV1 {gene}: trying {pdb_id}", flush=True)
            r = analyse(pdb_id, ref_seq, aligner)
            if r and (chosen is None or r["coverage"] > chosen["coverage"]):
                chosen = r
            if chosen and chosen["coverage"] > 0.85 and chosen["n_ligand_atoms"]:
                break
        if not chosen:
            print(f"  no usable structure for {gene}")
            continue
        summary[f"HSV1:{gene}"] = {k: chosen[k] for k in
                                   ("pdb_id", "chain", "coverage",
                                    "n_ligand_atoms", "n_contacts",
                                    "site_resolved")}
        summary[f"HSV1:{gene}"]["positions_with_features"] = len(chosen["features"])
        print(f"  chose {chosen['pdb_id']} chain {chosen['chain']}, "
              f"coverage {chosen['coverage']:.2f}, "
              f"{len(chosen['features'])} positions, "
              f"{chosen['n_ligand_atoms']} ligand atoms", flush=True)
        for pos, f in sorted(chosen["features"].items()):
            rows.append({
                "virus": "HSV1", "gene": gene, "position": pos,
                "rsa": round(f["rsa"], 4) if f["rsa"] is not None else "",
                "dist_to_ligand": round(f["dist_to_ligand"], 2)
                if f["dist_to_ligand"] is not None else "",
                "site_resolved": chosen["site_resolved"],
                "pdb_id": chosen["pdb_id"], "inherited": 0,
            })

        # HSV-2 inherits through a pairwise alignment to the HSV-1 reference.
        h2 = REF_DIR / f"HSV2_{gene}.faa"
        if not h2.exists():
            continue
        s2 = read_fasta_one(h2)
        aln = aligner.align(ref_seq, s2)[0]
        n_inh = 0
        for (rs, re_), (qs, _) in zip(*aln.aligned):
            for off in range(re_ - rs):
                f = chosen["features"].get(rs + off + 1)
                if not f:
                    continue
                rows.append({
                    "virus": "HSV2", "gene": gene, "position": qs + off + 1,
                    "rsa": round(f["rsa"], 4) if f["rsa"] is not None else "",
                    "dist_to_ligand": round(f["dist_to_ligand"], 2)
                    if f["dist_to_ligand"] is not None else "",
                    "site_resolved": chosen["site_resolved"],
                    "pdb_id": chosen["pdb_id"], "inherited": 1,
                })
                n_inh += 1
        summary[f"HSV2:{gene}"] = {"inherited_from": chosen["pdb_id"],
                                   "site_resolved": chosen["site_resolved"],
                                   "positions_with_features": n_inh}

    if not rows:
        print("no structural features produced")
        return 1
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "structure_features.tsv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    (REP / "structure_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(f"\n{len(rows)} structural feature rows written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
