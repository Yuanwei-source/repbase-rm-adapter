#!/usr/bin/env python3
"""
Convert RepBase EMBL ``*.ref`` records into a RepeatMasker ``-lib`` FASTA library.

Pipeline (per record):

  1. Parse ``ID``, ``KW`` tokens, ``OS`` species and the ``SQ`` sequence from EMBL.
  2. Classification: case-insensitive FIRST keyword that matches the Terrier
     ``repbase-to-repeatmasker.toml`` mapping table.  Unmatched records are kept
     as ``Unknown`` (listed separately, never force-guessed).
  3. Base normalization: ``x/X -> N``, ``u/U -> T``, ``o/O -> N``; uppercase;
     every other non-IUPAC character is dropped (there are none in practice).
     Every conversion event is logged.
  4. ID sanitization: ``: ( ) @`` are replaced with ``_``.
  5. De-duplication: identical sequence under the same ID is kept once; distinct
     sequences under the same ID are renamed ``ID_dup2``/``ID_dup3`` ...
  6. Emit a FASTA header ``>ID#class/subclass @species``.

Outputs (in addition to the FASTA library): ``id_mapping.tsv``,
``unmapped_records.tsv``, ``base_normalization.tsv`` and ``duplicate_report.tsv``,
written next to the output library (or to ``--out-dir`` when given).

The original RepBase files are never modified.
"""

from __future__ import annotations

import argparse
import collections
import glob
import os
import re
import sys
import tomllib

# Valid IUPAC nucleotide ambiguity codes (all are preserved as-is).
IUPAC = set("ACGTUNRYKMSWBDHV")


# --------------------------------------------------------------------------- #
# Mapping table
# --------------------------------------------------------------------------- #
def load_mapping(path: str) -> dict[str, str]:
    """Load the Terrier RepBase -> RepeatMasker TOML mapping.

    Keys are lower-cased for case-insensitive lookup.  Returns the lowercase
    key -> classification dict.  Raises if two keys collide case-insensitively
    onto different classifications.
    """
    with open(path, "rb") as f:
        raw = tomllib.load(f)

    mapping: dict[str, str] = {}
    for key, value in raw.items():
        if not isinstance(value, str):
            raise ValueError(f"mapping value for {key!r} is not a string")
        lkey = key.lower()
        if lkey in mapping and mapping[lkey] != value:
            raise ValueError(
                f"case-insensitive collision: {key!r} and an earlier key "
                f"map to {value!r} vs {mapping[lkey]!r}"
            )
        mapping[lkey] = value
    return mapping


# --------------------------------------------------------------------------- #
# EMBL parsing
# --------------------------------------------------------------------------- #
def parse_ref(path: str) -> list[dict]:
    """Parse one RepBase EMBL ``.ref`` file into a list of record dicts.

    Each record dict has ``id``, ``kw`` (list of semicolon-separated tokens),
    ``os`` (species) and ``seq`` (concatenated sequence).
    """
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    records = []
    for part in re.split(r"(?m)^//\s*$", text):
        m = re.search(r"(?m)^ID\s+(\S+)", part)
        if not m:
            continue
        rid = m.group(1)

        kw = []
        for line in re.findall(r"(?m)^KW\s+(.*)$", part):
            for tok in line.split(";"):
                tok = tok.strip().rstrip(".").strip()
                if tok:
                    kw.append(tok)

        os_m = re.search(r"(?m)^OS\s+(.+)$", part)
        os_ = os_m.group(1).strip() if os_m else None

        seq = ""
        sq_m = re.search(r"(?m)^SQ\s+Sequence\s+(\d+)\s+BP", part)
        if sq_m:
            for line in part[sq_m.end():].splitlines():
                line = line.strip()
                if not line:
                    continue
                m3 = re.match(r"^([A-Za-z\s]+?)\s*\d*\s*$", line)
                if m3:
                    seq += re.sub(r"\s+", "", m3.group(1))

        records.append(dict(id=rid, kw=kw, os=os_, seq=seq))
    return records


# --------------------------------------------------------------------------- #
# Normalization / classification / ID sanitization
# --------------------------------------------------------------------------- #
def normalize_bases(seq: str) -> tuple[str, collections.Counter]:
    """Normalize a sequence to uppercase IUPAC, mapping ``x/X->N``, ``u/U->T``,
    ``o/O->N`` and dropping any remaining non-IUPAC character.

    Returns ``(normalized_sequence, conversion_counter)``.
    """
    out: list[str] = []
    changes: collections.Counter = collections.Counter()
    for ch in seq:
        u = ch.upper()
        if u == "X":
            out.append("N")
            changes["x->N"] += 1
        elif u == "U":
            out.append("T")
            changes["u->T"] += 1
        elif u == "O":
            out.append("N")
            changes["o->N"] += 1
        elif u in IUPAC:
            out.append(u)
        else:
            changes[f"{ch!r}->DROP"] += 1
    return "".join(out), changes


def classify(kw: list[str], mapping: dict[str, str]) -> tuple[str, str | None]:
    """Return the RepeatMasker classification for a ``KW`` token list.

    The FIRST token that matches the mapping (case-insensitively) wins.  Returns
    ``(classification, matched_token)``; ``("Unknown", None)`` when nothing matches.
    """
    for tok in kw:
        val = mapping.get(tok.lower())
        if val is not None:
            return val, tok
    return "Unknown", None


def sanitize_id(rid: str) -> str:
    """Replace characters illegal in RepeatMasker FASTA IDs with ``_``."""
    return re.sub(r"[():@]", "_", rid)


def species_str(os_: str | None) -> str:
    """Normalize a species name for the FASTA header (spaces become ``_``)."""
    if not os_:
        return "Unknown"
    return os_.replace(" ", "_")


# --------------------------------------------------------------------------- #
# Main conversion
# --------------------------------------------------------------------------- #
def convert(input_dir: str, mapping_path: str, output_path: str,
            out_dir: str | None = None) -> dict:
    """Convert all ``*.ref`` files under ``input_dir`` into a RepeatMasker
    ``-lib`` FASTA library and audit reports.

    Returns a summary dict with the counts needed for validation.
    """
    mapping = load_mapping(mapping_path)
    files = sorted(glob.glob(os.path.join(input_dir, "*.ref")))
    if not files:
        raise FileNotFoundError(f"no *.ref files found under {input_dir}")

    all_recs: list[dict] = []
    for fp in files:
        all_recs.extend(parse_ref(fp))

    # Group records by original ID to detect duplicates.
    by_id: collections.OrderedDict = collections.OrderedDict()
    for r in all_recs:
        by_id.setdefault(r["id"], []).append(r)

    final_records: list[tuple[str, str, str, str, str]] = []  # new_id, orig_id, class, species, seq
    id_map_rows: list[tuple[str, str, str, str, str]] = []
    dup_report: list[tuple[str, str, int, str]] = []
    base_log: list[tuple[str, str, int]] = []
    unmapped: list[tuple[str, str]] = []
    orig_id_seen: collections.Counter = collections.Counter()

    for orig_id, recs in by_id.items():
        normed = []
        for r in recs:
            seq, changes = normalize_bases(r["seq"])
            if changes:
                for conv, n in changes.items():
                    base_log.append((r["id"], conv, n))
            normed.append((r, seq))

        # Group identical (normalized) sequences together.
        seq_groups: collections.OrderedDict = collections.OrderedDict()
        for r, seq in normed:
            seq_groups.setdefault(seq, []).append(r)

        # Classification uses the first record's KW (same ID shares KW in practice).
        classval, _matched_tok = classify(normed[0][0]["kw"], mapping)
        if classval == "Unknown":
            unmapped.append((orig_id, "; ".join(normed[0][0]["kw"])))

        idx = 0
        for seq, rlist in seq_groups.items():
            rep = rlist[0]
            idx += 1
            new_id = sanitize_id(orig_id) if idx == 1 else f"{sanitize_id(orig_id)}_dup{idx}"
            sp = species_str(rep["os"])
            final_records.append((new_id, orig_id, classval, sp, seq))
            for r in rlist:
                orig_id_seen[r["id"]] += 1
            if len(rlist) > 1:
                dup_report.append((orig_id, "exact_duplicate", len(rlist),
                                   f"kept {new_id}, dropped {len(rlist) - 1}"))
            id_map_rows.append((orig_id, new_id, classval, sp, f"seqgroup_{idx}"))

    # Post-sanitization ID uniqueness check.
    final_ids = collections.Counter(r[0] for r in final_records)
    collisions = {i: c for i, c in final_ids.items() if c > 1}
    if collisions:
        print(f"WARN: post-sanitization ID collisions: {collisions}", file=sys.stderr)

    # ---- write FASTA ----
    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        for new_id, _orig_id, classval, sp, seq in final_records:
            f.write(f">{new_id}#{classval} @{sp}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")

    # ---- write reports ----
    out_dir = out_dir or os.path.dirname(os.path.abspath(output_path))
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(out_dir, "id_mapping.tsv"), "w") as f:
        f.write("orig_id\tnew_id\tclass\tspecies\tnote\n")
        for orig_id, new_id, classval, sp, note in id_map_rows:
            f.write(f"{orig_id}\t{new_id}\t{classval}\t{sp}\t{note}\n")

    with open(os.path.join(out_dir, "unmapped_records.tsv"), "w") as f:
        f.write("orig_id\tkw\n")
        for orig_id, kw in unmapped:
            f.write(f"{orig_id}\t{kw}\n")

    with open(os.path.join(out_dir, "base_normalization.tsv"), "w") as f:
        f.write("orig_id\tconversion\tcount\n")
        for orig_id, conv, n in base_log:
            f.write(f"{orig_id}\t{conv}\t{n}\n")

    with open(os.path.join(out_dir, "duplicate_report.tsv"), "w") as f:
        f.write("orig_id\tkind\tn_records\tnote\n")
        for row in dup_report:
            f.write("\t".join(map(str, row)) + "\n")

    # ---- build summary ----
    class_hist = collections.Counter(r[2].split("/")[0] for r in final_records)
    base_totals: collections.Counter = collections.Counter()
    for _orig_id, conv, n in base_log:
        base_totals[conv] += n

    return {
        "input_records": len(all_recs),
        "unique_orig_ids": len(by_id),
        "duplicated_orig_ids": sum(1 for v in by_id.values() if len(v) > 1),
        "output_records": len(final_records),
        "dropped": len(all_recs) - len(final_records),
        "class_hist": class_hist,
        "unknown_unique_ids": len(unmapped),
        "unknown_records": class_hist.get("Unknown", 0),
        "base_totals": base_totals,
        "output_path": os.path.abspath(output_path),
        "out_dir": os.path.abspath(out_dir),
    }


def _print_summary(s: dict) -> None:
    print("\n================ SUMMARY ================")
    print(f"input records:        {s['input_records']}")
    print(f"unique original IDs:  {s['unique_orig_ids']}")
    print(f"duplicated original IDs: {s['duplicated_orig_ids']}")
    print(f"output FASTA records: {s['output_records']}")
    print(f"dropped (exact duplicates): {s['dropped']}")
    print("\nclassification (RepeatMasker top-level):")
    for k, v in s["class_hist"].most_common():
        print(f"  {k}: {v}")
    print(f"  Unknown (unique-ID groups): {s['unknown_unique_ids']}")
    print(f"  Unknown (output records):   {s['unknown_records']}")
    print("\nbase normalization events:")
    for conv, n in s["base_totals"].most_common():
        print(f"  {conv}: {n}")
    print(f"\nlibrary: {s['output_path']}")
    print(f"reports: {s['out_dir']}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert RepBase EMBL *.ref records into a RepeatMasker -lib FASTA library."
    )
    parser.add_argument("--input-dir", required=True,
                        help="directory containing RepBase *.ref files")
    parser.add_argument("--mapping", required=True,
                        help="path to repbase-to-repeatmasker.toml")
    parser.add_argument("--output", required=True,
                        help="destination FASTA library path")
    parser.add_argument("--out-dir", default=None,
                        help="directory for audit TSV reports (default: output directory)")
    args = parser.parse_args(argv)

    summary = convert(args.input_dir, args.mapping, args.output, args.out_dir)
    _print_summary(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
