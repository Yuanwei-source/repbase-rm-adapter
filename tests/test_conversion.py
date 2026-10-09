#!/usr/bin/env python3
"""
Synthetic-data tests for the RepBase -> RepeatMasker adapter.

These tests use short, hand-written EMBL records (NOT copied from RepBase) to
check the conversion behaviour:

* case-insensitive ``KW`` classification via the mapping table
* ``Unknown`` records are retained, never force-guessed
* same ID + same sequence -> de-duplicated; same ID + different sequence -> renamed
* ``x/u/o`` base normalization
* output FASTA ID uniqueness and record/class count consistency

Run with::

    python3 tests/test_conversion.py

or::

    python3 -m unittest discover -s tests -v
"""

from __future__ import annotations

import collections
import os
import shutil
import sys
import tempfile
import unittest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

import repbase_to_rm as r2r  # noqa: E402

MAPPING = os.path.join(REPO_ROOT, "mappings", "repbase-to-repeatmasker.toml")


def make_ref(records):
    """Build a synthetic RepBase ``.ref`` file from (id, kw, os, seq) tuples."""
    lines = []
    for rid, kw, os_, seq in records:
        lines.append(f"ID   {rid}")
        lines.append(f"KW   {kw}")
        if os_:
            lines.append(f"OS   {os_}")
        lines.append(f"SQ   Sequence {len(seq)} BP;")
        for i in range(0, len(seq), 60):
            lines.append("     " + seq[i:i + 60].lower())
        lines.append("//")
    return "\n".join(lines) + "\n"


class TestClassification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mapping = r2r.load_mapping(MAPPING)

    def test_gypsy(self):
        self.assertEqual(r2r.classify(["Gypsy"], self.mapping), ("LTR/Gypsy", "Gypsy"))

    def test_case_insensitive(self):
        self.assertEqual(r2r.classify(["gypsy"], self.mapping), ("LTR/Gypsy", "gypsy"))

    def test_first_match_wins(self):
        self.assertEqual(
            r2r.classify(["Mariner/Tc1", "DNA transposon", "Transposable Element"], self.mapping),
            ("DNA/TcMar", "Mariner/Tc1"),
        )

    def test_helitron_is_rc_not_dna(self):
        self.assertEqual(
            r2r.classify(["Helitron", "DNA transposon"], self.mapping),
            ("RC/Helitron", "Helitron"),
        )

    def test_penelope_is_ple(self):
        self.assertEqual(r2r.classify(["Penelope"], self.mapping), ("PLE", "Penelope"))

    def test_unknown_preserved(self):
        self.assertEqual(r2r.classify(["Totally Novel Element"], self.mapping), ("Unknown", None))

    def test_empty_kw_unknown(self):
        self.assertEqual(r2r.classify([], self.mapping), ("Unknown", None))


class TestNormalization(unittest.TestCase):
    def test_xuo(self):
        seq, changes = r2r.normalize_bases("acgtxuo")
        self.assertEqual(seq, "ACGTNTN")
        self.assertEqual(changes, collections.Counter({"x->N": 1, "u->T": 1, "o->N": 1}))

    def test_iupac_preserved_uppercased(self):
        seq, changes = r2r.normalize_bases("acgtrykmswbdhvn")
        self.assertEqual(seq, "ACGTRYKMSWBDHVN")
        self.assertEqual(changes, collections.Counter())


class TestSanitizeId(unittest.TestCase):
    def test_sanitize(self):
        self.assertEqual(r2r.sanitize_id("A@1"), "A_1")
        self.assertEqual(r2r.sanitize_id("tRNA-Leu-TTA(m)"), "tRNA-Leu-TTA_m_")
        self.assertEqual(r2r.sanitize_id("Scaffold_113:7041-7503"), "Scaffold_113_7041-7503")


class TestConvert(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_full_conversion(self):
        records = [
            ("G1", "Gypsy", "Anopheles gambiae", "acgtacgt"),
            ("G1", "Gypsy", "Anopheles gambiae", "acgtacgt"),       # exact duplicate
            ("G1", "Gypsy", "Anopheles gambiae", "acgtacgtacgt"),   # same ID, diff seq
            ("UNK1", "Totally Novel Element", "Homo sapiens", "tttt"),
            ("HEL1", "Helitron; DNA transposon; Transposable Element", "Mus musculus", "gggg"),
            ("NORM1", "Gypsy", "Danio rerio", "acgtxuo"),
        ]
        ref_path = os.path.join(self.tmp, "syn.ref")
        with open(ref_path, "w") as f:
            f.write(make_ref(records))

        out_fa = os.path.join(self.tmp, "out.fa")
        summary = r2r.convert(self.tmp, MAPPING, out_fa, self.tmp)

        self.assertEqual(summary["input_records"], 6)
        self.assertEqual(summary["unique_orig_ids"], 4)
        self.assertEqual(summary["duplicated_orig_ids"], 1)
        self.assertEqual(summary["output_records"], 5)
        self.assertEqual(summary["dropped"], 1)
        self.assertEqual(dict(summary["class_hist"]),
                         {"LTR": 3, "RC": 1, "Unknown": 1})
        self.assertEqual(summary["unknown_unique_ids"], 1)
        self.assertEqual(summary["unknown_records"], 1)
        self.assertEqual(dict(summary["base_totals"]),
                         {"x->N": 1, "u->T": 1, "o->N": 1})

        # Read the output FASTA and check headers + uniqueness.
        headers = []
        seqs = {}
        with open(out_fa) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    headers.append(line)
                    name = line[1:].split("#")[0]
                    seqs[name] = ""
                    continue
                if headers:
                    name = headers[-1][1:].split("#")[0]
                    seqs[name] += line

        ids = [h[1:].split("#")[0] for h in headers]
        self.assertEqual(len(ids), 5)
        self.assertEqual(len(set(ids)), 5, "FASTA IDs must be unique")

        # Record total == class total consistency.
        class_sum = sum(summary["class_hist"].values())
        self.assertEqual(class_sum, summary["output_records"])

        # Specific header + normalization checks.
        self.assertIn(">HEL1#RC/Helitron @Mus_musculus", headers)
        self.assertIn(">G1#LTR/Gypsy @Anopheles_gambiae", headers)
        self.assertIn(">G1_dup2#LTR/Gypsy @Anopheles_gambiae", headers)
        self.assertIn(">UNK1#Unknown @Homo_sapiens", headers)
        self.assertEqual(seqs.get("NORM1"), "ACGTNTN")


if __name__ == "__main__":
    unittest.main(verbosity=2)
