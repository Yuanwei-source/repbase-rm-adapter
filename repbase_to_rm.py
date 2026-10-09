#!/usr/bin/env python3
"""
Convert RepBase EMBL ``*.ref`` records into a RepeatMasker ``-lib`` FASTA library.

Pipeline (per record):

  1. Parse ``ID``, ``KW`` tokens, ``OS`` species and the ``SQ`` sequence from EMBL.
     Each record is validated: the ``SQ`` field must exist, the sequence must be
     non-empty, and its parsed length must equal the declared ``N BP`` count
     (relaxable with ``--allow-length-mismatch``, which logs instead of aborting).
  2. Classification: case-insensitive FIRST keyword that matches the Terrier
     ``repbase-to-repeatmasker.toml`` mapping table.  Unmatched records are kept
     as ``Unknown`` (listed separately, never force-guessed).  Classification is
     done per distinct normalized sequence group, so records that share an ID but
     differ in sequence (and ``KW``) each get their own class.
  3. Base normalization: ``x/X -> N``, ``u/U -> T``, ``o/O -> N``; uppercase.
     Legal IUPAC ambiguity codes are preserved; every other character is an
     error, never a silent drop (a drop would shorten the sequence and shift all
     downstream coordinates).  Every conversion event is logged.
  4. ID sanitization: ``: ( ) @`` are replaced with ``_``.
  5. De-duplication: identical sequence under the same ID is kept once; distinct
     sequences under the same ID are renamed ``ID_dup2``/``ID_dup3`` ...
  6. Emit a FASTA header ``>ID#class/subclass @species``.

The converter refuses to write a library it cannot trust: malformed ``SQ``
blocks, illegal bases and post-sanitization FASTA ID collisions all abort with a
non-zero exit code.

When several records share an identical normalized sequence but carry different
classifications, one deterministic label is chosen (a concrete class beats
``Unknown``, a class with a subclass beats a bare class, and otherwise the first
record in input order wins) and every such event is listed in
``class_conflicts.tsv`` and counted in the summary -- nothing is overwritten
silently.

Outputs (in addition to the FASTA library): ``id_mapping.tsv``,
``unmapped_records.tsv``, ``base_normalization.tsv``, ``duplicate_report.tsv``
and ``class_conflicts.tsv`` (written only when conflicts occurred), written next
to the output library (or to ``--out-dir`` when given).

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


class ConversionError(Exception):
    """Raised when input data would produce an invalid or ambiguous library.

    The converter never writes a partially-validated library and never silently
    repairs data; the caller aborts with a non-zero exit code instead.
    """


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
def parse_ref(path: str,
              allow_length_mismatch: bool = False) -> tuple[list[dict], list[tuple]]:
    """Parse one RepBase EMBL ``.ref`` file into a list of record dicts.

    Each record dict has ``id``, ``kw`` (list of semicolon-separated tokens),
    ``os`` (species), ``src`` (source file) and ``seq`` (concatenated sequence).

    Returns ``(records, length_anomalies)``.  Raises :class:`ConversionError`
    when a record has no ``SQ`` field or an empty sequence, or when the parsed
    sequence length differs from the declared ``N BP`` count and
    ``allow_length_mismatch`` is false.
    """
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    records = []
    length_anomalies: list[tuple] = []
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

        sq_m = re.search(r"(?m)^SQ\s+Sequence\s+(\d+)\s+BP", part)
        if not sq_m:
            raise ConversionError(
                f"{path}: record {rid}: no 'SQ Sequence N BP' field found; "
                f"refusing to emit a record without a sequence"
            )

        seq = ""
        for line in part[sq_m.end():].splitlines():
            line = line.strip()
            if not line:
                continue
            m3 = re.match(r"^([A-Za-z\s]+?)\s*\d*\s*$", line)
            if m3:
                seq += re.sub(r"\s+", "", m3.group(1))

        if not seq:
            raise ConversionError(
                f"{path}: record {rid}: the SQ field declares a sequence but no "
                f"nucleotide characters could be parsed"
            )

        declared = int(sq_m.group(1))
        if len(seq) != declared:
            if not allow_length_mismatch:
                raise ConversionError(
                    f"{path}: record {rid}: SQ declares {declared} BP but "
                    f"{len(seq)} characters were parsed; refusing to emit an "
                    f"inconsistent record (pass --allow-length-mismatch to accept "
                    f"and log these instead)"
                )
            length_anomalies.append((path, rid, declared, len(seq)))

        records.append(dict(id=rid, kw=kw, os=os_, src=path, seq=seq))
    return records, length_anomalies


# --------------------------------------------------------------------------- #
# Normalization / classification / ID sanitization
# --------------------------------------------------------------------------- #
def normalize_bases(seq: str, rid: str | None = None,
                    src: str | None = None) -> tuple[str, collections.Counter]:
    """Normalize a sequence to uppercase IUPAC, mapping ``x/X->N``, ``u/U->T``
    and ``o/O->N``.

    Legal IUPAC ambiguity codes are preserved.  Any other character raises
    :class:`ConversionError`: dropping it would silently shorten the sequence
    and shift every downstream coordinate.

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
            where = f"{src or '<input>'}: record {rid}: " if rid is not None else ""
            raise ConversionError(
                f"{where}illegal nucleotide {ch!r} in sequence; only IUPAC bases "
                f"(plus x/u/o variants) are accepted, and illegal characters are "
                f"never silently dropped"
            )
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


def resolve_class_conflict(candidates):
    """Choose one classification when identical sequences disagree.

    ``candidates`` is an iterable of ``(position, classvalue)`` pairs in input
    order.  The rule is deterministic and documented:

    1. a concrete classification beats ``Unknown``;
    2. a class with a subclass (``LTR/Gypsy``) beats a bare class (``LTR``);
    3. otherwise the first record in input order wins.

    Callers are expected to record every conflict; this helper only picks.
    """
    def rank(item):
        position, classval = item
        return (1 if classval == "Unknown" else 0,
                0 if "/" in classval else 1,
                position)

    return min(candidates, key=rank)[1]


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
            out_dir: str | None = None,
            allow_length_mismatch: bool = False) -> dict:
    """Convert all ``*.ref`` files under ``input_dir`` into a RepeatMasker
    ``-lib`` FASTA library and audit reports.

    Returns a summary dict with the counts needed for validation.  Raises
    :class:`ConversionError` instead of writing an invalid or ambiguous library.
    """
    mapping = load_mapping(mapping_path)
    files = sorted(glob.glob(os.path.join(input_dir, "*.ref")))
    if not files:
        raise FileNotFoundError(f"no *.ref files found under {input_dir}")

    all_recs: list[dict] = []
    length_anomalies: list[tuple] = []
    for fp in files:
        recs, anomalies = parse_ref(fp, allow_length_mismatch=allow_length_mismatch)
        all_recs.extend(recs)
        length_anomalies.extend(anomalies)

    # Group records by original ID to detect duplicates.
    by_id: collections.OrderedDict = collections.OrderedDict()
    for r in all_recs:
        by_id.setdefault(r["id"], []).append(r)

    final_records: list[tuple[str, str, str, str, str]] = []  # new_id, orig_id, class, species, seq
    id_map_rows: list[tuple[str, str, str, str, str]] = []
    dup_report: list[tuple[str, str, int, str]] = []
    base_log: list[tuple[str, str, int]] = []
    unmapped: list[tuple[str, str, str]] = []
    conflict_rows: list[tuple] = []
    orig_id_seen: collections.Counter = collections.Counter()

    for orig_id, recs in by_id.items():
        normed = []
        for r in recs:
            seq, changes = normalize_bases(r["seq"], rid=r["id"], src=r["src"])
            if changes:
                for conv, n in changes.items():
                    base_log.append((r["id"], conv, n))
            normed.append((r, seq))

        # Group identical (normalized) sequences together.
        seq_groups: collections.OrderedDict = collections.OrderedDict()
        for r, seq in normed:
            seq_groups.setdefault(seq, []).append(r)

        idx = 0
        for seq, rlist in seq_groups.items():
            rep = rlist[0]
            # Classify every distinct normalized sequence independently from its
            # own KW: same-ID records that differ in sequence may legitimately
            # belong to different classes.
            classed = [(pos, classify(r["kw"], mapping)[0])
                       for pos, r in enumerate(rlist)]
            group_classes = {c for _pos, c in classed}
            if len(group_classes) == 1:
                classval = classed[0][1]
            else:
                # Records sharing an identical sequence disagree on the class.
                # Resolve deterministically and record every occurrence; never
                # overwrite one annotation with another silently.
                classval = resolve_class_conflict(classed)
                counts = collections.Counter(c for _pos, c in classed)
                conflict_rows.append((
                    orig_id,
                    len(seq),
                    classval,
                    " | ".join(f"{c}:{n}" for c, n in sorted(counts.items())),
                    ", ".join(sorted({os.path.basename(r["src"]) for r in rlist})),
                ))

            idx += 1
            new_id = sanitize_id(orig_id) if idx == 1 else f"{sanitize_id(orig_id)}_dup{idx}"
            sp = species_str(rep["os"])
            if classval == "Unknown":
                unmapped.append((orig_id, new_id, "; ".join(rep["kw"])))
            final_records.append((new_id, orig_id, classval, sp, seq))
            for r in rlist:
                orig_id_seen[r["id"]] += 1
            if len(rlist) > 1:
                dup_report.append((orig_id, "exact_duplicate", len(rlist),
                                   f"kept {new_id}, dropped {len(rlist) - 1}"))
            id_map_rows.append((orig_id, new_id, classval, sp, f"seqgroup_{idx}"))

    # Post-sanitization ID uniqueness check.  Sanitization and the ``_dupN``
    # renaming can collide with another entry's natural ID, which would produce a
    # library with duplicate FASTA names.
    final_ids = collections.Counter(r[0] for r in final_records)
    collisions = {i: c for i, c in final_ids.items() if c > 1}
    if collisions:
        detail = ", ".join(f"{i!r} x{c}" for i, c in sorted(collisions.items()))
        raise ConversionError(
            f"post-sanitization FASTA ID collisions ({detail}); different input "
            f"entries would map onto the same output name. Refusing to write a "
            f"library with duplicate IDs."
        )

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
        f.write("orig_id\tnew_id\tkw\n")
        for orig_id, new_id, kw in unmapped:
            f.write(f"{orig_id}\t{new_id}\t{kw}\n")

    with open(os.path.join(out_dir, "base_normalization.tsv"), "w") as f:
        f.write("orig_id\tconversion\tcount\n")
        for orig_id, conv, n in base_log:
            f.write(f"{orig_id}\t{conv}\t{n}\n")

    with open(os.path.join(out_dir, "duplicate_report.tsv"), "w") as f:
        f.write("orig_id\tkind\tn_records\tnote\n")
        for row in dup_report:
            f.write("\t".join(map(str, row)) + "\n")

    if conflict_rows:
        with open(os.path.join(out_dir, "class_conflicts.tsv"), "w") as f:
            f.write("orig_id\tseq_len\tchosen_class\tcandidates\tsources\n")
            for row in conflict_rows:
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
        "unknown_unique_ids": len({o for o, _n, _k in unmapped}),
        "unknown_records": class_hist.get("Unknown", 0),
        "base_totals": base_totals,
        "length_anomalies": length_anomalies,
        "class_conflicts": conflict_rows,
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
    if s.get("length_anomalies"):
        print(f"\nSQ length anomalies accepted (--allow-length-mismatch): "
              f"{len(s['length_anomalies'])}")
        for path, rid, declared, actual in s["length_anomalies"]:
            print(f"  {os.path.basename(path)}: {rid}: declared {declared}, parsed {actual}")
    if s.get("class_conflicts"):
        print(f"\nidentical-sequence class conflicts resolved: "
              f"{len(s['class_conflicts'])} (see class_conflicts.tsv)")
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
    parser.add_argument("--allow-length-mismatch", action="store_true",
                        help="accept records whose parsed sequence length differs "
                             "from the SQ-declared BP count (logged in the summary "
                             "instead of aborting)")
    args = parser.parse_args(argv)

    try:
        summary = convert(args.input_dir, args.mapping, args.output, args.out_dir,
                          allow_length_mismatch=args.allow_length_mismatch)
    except ConversionError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    _print_summary(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
