"""Guard for the gene-name matcher. Run directly: python scripts/test_match_gene.py

The matcher decides which GenBank CDS records enter the analysis, so a change
to it silently changes the isolate set and every downstream number. Two
failure directions matter and they pull against each other:

  Too strict -- a record that plainly is UL23 is skipped. This cost 244
  accessions before the bare-locus-name pass was added, and hit HSV-2 hardest.

  Too loose  -- UL5 captures UL50 to UL56, or a helicase-primase subunit of
  unknown identity is guessed at. Those six UL5x genes appear in nearly every
  record, so a substring match would contaminate the alignment with unrelated
  proteins and the error would not be obvious in any summary statistic.

Both directions are covered below.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p1_02_gate2 import match_gene                             # noqa: E402

CASES = [
    # --- descriptive product names
    ({"protein": "thymidine kinase"}, "UL23"),
    ({"protein": "DNA polymerase"}, "UL30"),
    ({"protein": "helicase-primase helicase subunit"}, "UL5"),
    ({"protein": "helicase-primase primase subunit"}, "UL52"),
    # --- explicit gene qualifiers and their aliases
    ({"gene": "ul23"}, "UL23"),
    ({"gene": "tk"}, "UL23"),
    ({"gene": "ul30"}, "UL30"),
    ({"gene": "pol"}, "UL30"),
    ({"gene": "ul5"}, "UL5"),
    ({"gene": "ul52"}, "UL52"),
    # --- bare locus name as the product, no gene qualifier
    ({"protein": "UL23"}, "UL23"),
    ({"protein": "UL30"}, "UL30"),
    ({"protein": "UL5"}, "UL5"),
    ({"protein": "UL52"}, "UL52"),
    ({"protein": "ul23"}, "UL23"),
    # --- the six neighbours a substring match would swallow
    ({"protein": "UL50"}, None),
    ({"protein": "UL51"}, None),
    ({"protein": "UL53"}, None),
    ({"protein": "UL54"}, None),
    ({"protein": "UL55"}, None),
    ({"protein": "UL56"}, None),
    ({"protein": "tegument protein UL51"}, None),
    # --- genes of the same complex or family that must not be captured
    ({"protein": "UL8"}, None),
    ({"protein": "UL9"}, None),
    ({"protein": "helicase-primase subunit"}, None),
    ({"protein": "DNA helicase/primase complex-associated protein"}, None),
    ({"protein": "DNA polymerase processivity subunit"}, None),
    ({"protein": "DNA polymerase accessory protein"}, None),
    # --- an explicit gene qualifier is authoritative and must not fall through
    #     to a protein-name match that would steal a paralog
    ({"gene": "ul9", "protein": "dna helicase"}, None),
    ({"gene": "ul8", "protein": "thymidine kinase"}, None),
    # --- nothing to go on
    ({}, None),
    ({"protein": ""}, None),
]


def main():
    failures = []
    for q, want in CASES:
        got = match_gene(q)
        if got != want:
            failures.append(f"  {q!r} -> {got!r}, expected {want!r}")
    print(f"{len(CASES)} cases")
    if failures:
        print(f"\n{len(failures)} FAILURE(S):", file=sys.stderr)
        for f in failures:
            print(f, file=sys.stderr)
        return 1
    print("all cases pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
