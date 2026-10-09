#!/usr/bin/env python3
"""
Synthetic-data tests for the RepBase -> RepeatMasker adapter.

These tests use short, hand-written EMBL records (NOT copied from RepBase) to
check the conversion behaviour:

* case-insensitive ``KW`` classification via the mapping table
* ``Unknown`` records are retained, never force-guessed
* same ID + same sequence -> de-duplicated; same ID + different sequence -> renamed
* same ID + different sequence + different ``KW`` -> classified independently
* same ID + same sequence + conflicting ``KW`` -> resolved deterministically
  (concrete beats ``Unknown``, subclass beats bare class, else input order) and
  every conflict is written to ``class_conflicts.tsv``
* ``x/u/o`` base normalization
* output FASTA ID uniqueness and record/class count consistency
* malformed input (missing ``SQ``, empty sequence, declared-length mismatch,
  illegal bases) and post-sanitization ID collisions abort instead of writing
  a plausible-but-wrong library
* cross-validation reports no-shared-names instead of dividing by zero

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


class TestInputValidation(unittest.TestCase):
    """Malformed input must abort instead of producing a plausible library."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def _write(self, text, name="syn.ref"):
        p = os.path.join(self.tmp, name)
        with open(p, "w") as f:
            f.write(text)
        return p

    def _convert(self, **kw):
        return r2r.convert(self.tmp, MAPPING, os.path.join(self.tmp, "out.fa"),
                           self.tmp, **kw)

    def test_post_sanitization_collision_aborts(self):
        # 'A@1' and 'A_1' both sanitize to the same FASTA ID.
        self._write(make_ref([("A@1", "Gypsy", "X", "acgt"),
                              ("A_1", "Gypsy", "X", "tttt")]))
        with self.assertRaises(r2r.ConversionError):
            self._convert()

    def test_generated_dup_suffix_collision_aborts(self):
        # ID 'A' with two distinct sequences emits 'A' and 'A_dup2', which
        # collides with the natural entry 'A_dup2'.
        self._write(make_ref([("A", "Gypsy", "X", "acgt"),
                              ("A", "Gypsy", "X", "tttt"),
                              ("A_dup2", "Gypsy", "X", "cccc")]))
        with self.assertRaises(r2r.ConversionError):
            self._convert()

    def test_missing_sq_aborts(self):
        self._write("ID   NOSQ\nKW   Gypsy\nOS   X\n//\n")
        with self.assertRaises(r2r.ConversionError) as cm:
            self._convert()
        self.assertIn("NOSQ", str(cm.exception))

    def test_empty_sequence_aborts(self):
        self._write("ID   EMPTY\nKW   Gypsy\nOS   X\nSQ   Sequence 10 BP;\n//\n")
        with self.assertRaises(r2r.ConversionError) as cm:
            self._convert()
        self.assertIn("EMPTY", str(cm.exception))

    def test_declared_length_mismatch_aborts_by_default(self):
        self._write("ID   SHORT\nKW   Gypsy\nOS   X\nSQ   Sequence 999 BP;\n     acgt\n//\n")
        with self.assertRaises(r2r.ConversionError) as cm:
            self._convert()
        self.assertIn("SHORT", str(cm.exception))

    def test_declared_length_mismatch_allowed_is_logged(self):
        self._write("ID   SHORT\nKW   Gypsy\nOS   X\nSQ   Sequence 999 BP;\n     acgt\n//\n")
        summary = self._convert(allow_length_mismatch=True)
        self.assertEqual(len(summary["length_anomalies"]), 1)
        self.assertEqual(summary["length_anomalies"][0][1], "SHORT")

    def test_illegal_base_aborts(self):
        self._write(make_ref([("BADB", "Gypsy", "X", "acgtz")]))
        with self.assertRaises(r2r.ConversionError) as cm:
            self._convert()
        self.assertIn("'z'", str(cm.exception))

    def test_illegal_character_is_never_dropped(self):
        with self.assertRaises(r2r.ConversionError):
            r2r.normalize_bases("ACGTQ")


class TestPerSequenceGroupClassification(unittest.TestCase):
    """Same ID + different sequence must be classified from its own KW."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def _run(self, records):
        with open(os.path.join(self.tmp, "syn.ref"), "w") as f:
            f.write(make_ref(records))
        out_fa = os.path.join(self.tmp, "out.fa")
        summary = r2r.convert(self.tmp, MAPPING, out_fa, self.tmp)
        with open(out_fa) as f:
            headers = [line.strip() for line in f if line.startswith(">")]
        return summary, headers

    def test_distinct_sequences_keep_their_own_classification(self):
        summary, headers = self._run([
            ("M1", "Gypsy", "X", "aaaa"),
            ("M1", "Helitron; DNA transposon", "X", "tttt"),
        ])
        self.assertIn(">M1#LTR/Gypsy @X", headers)
        self.assertIn(">M1_dup2#RC/Helitron @X", headers)
        self.assertEqual(dict(summary["class_hist"]), {"LTR": 1, "RC": 1})

    def test_conflicting_class_for_identical_sequence_is_resolved_and_recorded(self):
        # Deterministic rule: both classes carry a subclass, so the first record
        # in input order wins.  The conflict is recorded, never silently dropped.
        summary, headers = self._run([("C1", "Gypsy", "X", "aaaa"),
                                      ("C1", "Helitron", "X", "aaaa")])
        self.assertEqual(headers, [">C1#LTR/Gypsy @X"])
        self.assertEqual(len(summary["class_conflicts"]), 1)
        self.assertEqual(summary["class_conflicts"][0][2], "LTR/Gypsy")
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "class_conflicts.tsv")))

    def test_unknown_loses_to_a_concrete_classification(self):
        summary, headers = self._run([("C2", "NoSuchKeyword", "X", "aaaa"),
                                      ("C2", "Gypsy", "X", "aaaa")])
        self.assertEqual(headers, [">C2#LTR/Gypsy @X"])
        self.assertEqual(summary["class_conflicts"][0][2], "LTR/Gypsy")

    def test_bare_class_loses_to_a_subclass(self):
        summary, headers = self._run([("C3", "DNA transposon", "X", "aaaa"),
                                      ("C3", "Mariner/Tc1", "X", "aaaa")])
        self.assertEqual(headers, [">C3#DNA/TcMar @X"])
        self.assertEqual(summary["class_conflicts"][0][2], "DNA/TcMar")


class TestCrossValidate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import importlib.util
        path = os.path.join(REPO_ROOT, "scripts", "cross_validate.py")
        spec = importlib.util.spec_from_file_location("cross_validate", path)
        cls.cv = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.cv)

    def test_no_shared_names_returns_error(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        a, b = os.path.join(d, "a.lib"), os.path.join(d, "b.lib")
        with open(a, "w") as f:
            f.write(">X1#LTR/Gypsy @A\nacgt\n")
        with open(b, "w") as f:
            f.write(">Y1#LTR/Gypsy @A\nacgt\n")
        self.assertNotEqual(self.cv.main(["--old-lib", a, "--new-lib", b]), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
