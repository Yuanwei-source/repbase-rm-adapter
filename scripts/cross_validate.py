#!/usr/bin/env python3
"""
Cross-validate a converted RepeatMasker library against a reference library.

Both libraries use the RepeatMasker FASTA header convention::

    >name#class/subclass @species

This script reports, for the names shared between the two libraries, the
top-level class agreement and the full class/subclass agreement.  It is useful
for checking that a converted RepBase library reproduces the classifications of
the official RepeatMasker Edition (or any other reference library).

Usage::

    python cross_validate.py --old-lib ref.lib --new-lib converted.fa
"""

from __future__ import annotations

import argparse
import collections
import sys


def parse_lib(path: str) -> dict[str, str]:
    """Return ``{sequence_id: classification}`` for a RepeatMasker-style FASTA."""
    d: dict[str, str] = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                hdr = line[1:].strip()
                name = hdr.split("#")[0]
                rest = hdr.split("#", 1)[1] if "#" in hdr else ""
                cls_part = rest.split("@")[0].strip() if rest else ""
                cls = cls_part.split()[0] if cls_part else ""
                d[name] = cls
    return d


def top_level(cls: str) -> str:
    return cls.split("/")[0] if cls else ""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Cross-validate a converted library against a reference library."
    )
    parser.add_argument("--old-lib", required=True, help="reference library (e.g. 2018 RM Edition)")
    parser.add_argument("--new-lib", required=True, help="converted library to validate")
    args = parser.parse_args(argv)

    old = parse_lib(args.old_lib)
    new = parse_lib(args.new_lib)
    print(f"reference library: {len(old)} entries")
    print(f"converted library: {len(new)} entries")

    shared = set(old) & set(new)
    print(f"shared names:      {len(shared)}")

    agree_top = 0
    agree_full = 0
    disagree: list[tuple[str, str, str]] = []
    for name in sorted(shared):
        oc, nc = old[name], new[name]
        if top_level(oc) == top_level(nc):
            agree_top += 1
        if oc == nc:
            agree_full += 1
        else:
            disagree.append((name, oc, nc))

    print(f"\ntop-level class agreement:   {agree_top}/{len(shared)} "
          f"({100 * agree_top / len(shared):.2f}%)")
    print(f"full class/subclass agreement: {agree_full}/{len(shared)} "
          f"({100 * agree_full / len(shared):.2f}%)")
    print(f"disagreements: {len(disagree)}")

    cat: collections.Counter = collections.Counter()
    for name, oc, nc in disagree:
        if nc == "Unknown":
            cat["new=Unknown"] += 1
        elif top_level(oc) != top_level(nc):
            cat[f"top-differ {oc}->{nc}"] += 1
        else:
            cat[f"subclass-differ {oc}->{nc}"] += 1
    print("\ndisagreement categories:")
    for k, v in cat.most_common(40):
        print(f"  {k}: {v}")
    print("\nsample disagreements (name: reference -> converted):")
    for name, oc, nc in disagree[:50]:
        print(f"  {name}: {oc} -> {nc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
