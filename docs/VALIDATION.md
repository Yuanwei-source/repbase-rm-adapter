# Validation Report

**English** · [简体中文](VALIDATION.zh-CN.md)

This document records how the RepBase-to-RepeatMasker adapter was validated, the
commands used, the results, and the boundaries of what was and was not tested.

**Nature of this validation: functional compatibility validation.** It
establishes that the converter produces a well-formed, unambiguous, complete
RepeatMasker library and that RepeatMasker can consume it correctly. It is **not**
a genome-wide accuracy benchmark, and it does **not** establish that this library
is superior to the curated 2018 RepeatMasker Edition for any particular genome.

---

## 1. Validation scope

| Item | Value |
| --- | --- |
| Converter | `repbase_to_rm.py` (this repository) |
| Input | RepBase 31.09 EMBL, 39 main `.ref` files (appendix excluded) |
| Input location | local; not distributed |
| Mapping | `mappings/repbase-to-repeatmasker.toml` (Terrier, unmodified) |
| Python | 3.13.13 (standard library only) |
| RepeatMasker | 4.2.4 |
| RMBlast | 2.17.1 |
| BLAST+ | `makeblastdb` (`-parse_seqids`) |
| Date | 2026-10-09 |
| Reference for cross-check | RepBase 2018 RepeatMasker Edition (`repbase_2018.lib`, 49,067 entries) |

**In scope:** record parsing, classification mapping, identifier handling,
base normalization, output format, RepeatMasker compatibility, and sequence
integrity against the official RepBase 31.09 FASTA release.

**Out of scope:** genome-wide sensitivity/specificity of the resulting
annotation; taxonomic restriction of the library; comparison against
sequence-based classifiers; the `appendix/` records (deliberately excluded
because they share 7,441 identifiers with the main files).

The validation is organised in three independent parts:

| Part | Sections | Evidence |
| --- | --- | --- |
| Conversion integrity | §2–§4 | record accounting, classification distribution, duplicate and base handling, official-FASTA sequence reconciliation, CLI reproducibility |
| RepeatMasker compatibility | §5 | `makeblastdb` build and `-lib` functional test with positive controls |
| Real-genome evaluation | §6 | three-library comparison on a real 10 Mb window |

---

## 2. Conversion summary

Command:

```bash
python repbase_to_rm.py \
    --input-dir /path/to/RepBase31.09.embl \
    --mapping   mappings/repbase-to-repeatmasker.toml \
    --output    repbase_31.09.rm.fa \
    --out-dir   reports \
    --allow-length-mismatch
```

`--allow-length-mismatch` is required for RepBase 31.09: five records in
`humsub.ref` (`AluY`, `AluYa5`, `AluYa8`, `AluYb8`, `AluYb9`) declare an `SQ`
base-pair count that disagrees with the length of their sequence lines. Without
the flag the converter aborts on such a record instead of silently writing a
truncated or padded sequence.

### 2.1 Record accounting

| Item | Count |
| --- | --- |
| Input EMBL records | 126,885 |
| Unique original IDs | 124,866 |
| Duplicated original IDs (appearing >1×) | 1,556 (3,575 records) |
| — of which fully identical (all records same sequence) | 1,513 (3,441 records) → drop 1,928 |
| — of which mixed (same ID, differing sequences) | 43 (134 records) → keep 107 distinct sequences, drop 27 internal exact duplicates |
| Renamed records added (`_dup2` … `_dup13`) | 64 |
| Output FASTA records | **124,930** |

Reconciliation identities:

```
126,885 − 1,955 (exact duplicates dropped) = 124,930
124,866 (unique IDs) + 64 (renamed)         = 124,930
```

`_dup2` … `_dup13`: the identifier with the most distinct sequences is
`Contig_26` (13 sequences).

### 2.2 Classification distribution (record-level, directly from the final FASTA)

Counts below are tallied directly from the FASTA headers of the delivered
library, so the table sums exactly to the record total.

| Top-level class | Records |
| --- | --- |
| LTR | 75,260 |
| DNA | 27,740 |
| LINE | 13,512 |
| SINE | 2,930 |
| RC | 2,198 |
| PLE | 1,534 |
| Unknown | 847 |
| Satellite | 759 |
| Structural_RNA | 139 |
| Other | 11 |
| **Total** | **124,930** |

### 2.3 Classification coverage

| Metric | Value |
| --- | --- |
| Classified records | 124,083 |
| Unclassified records | 847 |
| Coverage (records) | **99.32%** |
| Distinct original entries that are unclassified | 823 |

The 847 unclassified **records** correspond to 823 distinct unclassified
**entries**: a few of those entries are represented by more than one distinct
sequence and therefore also by a renamed `_dupN` record carrying the same
`Unknown` class, which accounts for the remaining 24 records. Both figures are
reported for transparency; the record-level figure (99.32%) is the one
consistent with the 124,930-record output.

### 2.4 Composition of the `Unknown` records

Of the 823 distinct unclassified entries:

- **Empty `KW` field** — 93 entries (117 output records): no classification
  information at all.
- **Generic descriptors, family undeterminable** — 548 entries whose `KW`
  contains one of the generic phrases `Transposable Element` (231 entries),
  `Nonautonomous` (254), `Repetitive sequence` (121), `conserved` (117),
  `Repetitive element` (106), `Multicopy gene` (26), `Repeat region` (11). These
  counts are entries containing the phrase and overlap freely; the 548 figure is
  the number of distinct entries matching at least one of them. Such entries are
  honestly marked `Unknown`, never force-guessed.
- **Named but unmapped** — the remaining 182 entries, whose `KW` names an element
  or family that the Terrier mapping does not cover (`SINE element`; structural
  RNA variants such as `Transfer RNA` / `Small nuclear RNA`; LTR families such
  as `Loki-*-LTR` / `CoeEFV-LTR`; `MITE`; `FLA`; `Telomeric repeat`; and entries
  whose `KW` is only an element name such as `IKIRARA1` or `MINIME_DN`). A
  recorded synonym supplement could recover part of these (raising coverage
  towards ≈99.6%), but per scope they are kept as `Unknown` rather than silently
  mapped.

The full list is in `unmapped_records.tsv` (823 distinct IDs, 847 records).

### 2.5 Base normalization (recorded, not blanket-replaced)

| Conversion | Characters | Records | Note |
| --- | --- | --- | --- |
| `x/X` → `N` | 37,502 | 221 | RepBase masked/unknown positions |
| `u/U` → `T` | 8 | 4 | uracil residues in DNA consensus |
| `o/O` → `N` | 1 | 1 | single illegal character in `Merlin-5N2_OT` |

All other non-ACGT characters are legal IUPAC ambiguity codes
(`N/R/Y/K/M/S/W/B/D/H/V`, ≈430,000 characters) and are preserved unchanged.
Every event is logged in `base_normalization.tsv`.

### 2.6 Identifier sanitization

166 records contained characters illegal in RepeatMasker FASTA identifiers
(`:`, `(`, `)`, `@`); these are replaced with `_` (e.g. `A@1` → `A_1`,
`tRNA-Leu-TTA(m)` → `tRNA-Leu-TTA_m_`). After sanitization, the output contains
**zero duplicate identifiers**; a collision is a hard error, so the converter
stops instead of writing duplicate FASTA identifiers.

### 2.7 Input validation and classification-conflict resolution

The converter validates its input before writing anything, and resolves the
ambiguity present in the real release by fixed rules rather than by read order.

**Hard checks — a failure aborts the run and no library is written:**

| Check | RepBase 31.09 |
| --- | --- |
| Record has an `SQ` line | 0 missing |
| Sequence is non-empty | 0 empty |
| Declared `SQ` length matches the sequence lines | 5 disagreements (below) |
| Sequence contains only IUPAC characters | 0 illegal characters |
| No duplicate identifier after sanitization | 0 collisions |

The five declared-length disagreements are properties of the released data, not
of the parser: the `SQ` line and the sequence lines disagree with each other
(`humsub.ref`: `AluY` 281/282, `AluYa5` 281/282, `AluYa8` 280/281, `AluYb8`
288/289, `AluYb9` 282/289). Because aborting would make the real release
unconvertible, `--allow-length-mismatch` accepts them explicitly and prints every
accepted anomaly in the run summary; without the flag the run stops rather than
silently writing a truncated or padded sequence. A missing or empty sequence, an
illegal base, and an identifier collision are always hard errors: the converter
stops instead of writing a plausible-but-wrong library.

**Classification conflicts.** When one identifier covers several records whose
normalized sequences are identical but whose `KW` fields map to different
classes, the class is resolved deterministically: a concrete class beats
`Unknown`, a class with a subclass beats a bare top-level class, and otherwise
the first record in input order wins. All 62 conflicts in RepBase 31.09 are
written to `class_conflicts.tsv` with their candidate classes and source files;
none is resolved silently. When one identifier covers records with *different*
sequences, each sequence is classified from the `KW` of its own record instead
of inheriting another record's class, and `id_mapping.tsv` records which
sequence group each output record came from (`seqgroup_N`).

**Interpretation limit.** Where one normalized sequence has several candidate
classes, the program resolves the conflict with a deterministic heuristic. That
rule makes the result reproducible; it is **not** an independent biological
verification of classification correctness. Every conflict is recorded in
`class_conflicts.tsv` for further review by the user.

**How far the conflicts reach.** Of the 62 conflicts, **35** are between
candidates in the same top-level class (e.g. `SINE` vs `SINE/tRNA`, `DNA` vs
`DNA/MULE`); **27** span different top-level labels, but **26** of those pair a
concrete class with `Unknown` rather than disagreeing about biology. Only one
record, `TS2`, is a genuine disagreement between two distinct concrete top-level
classes (`SINE` in `dcotrep.ref` vs `LINE/RTE` in `edcotrep.ref`, the latter last
updated in 2021); no other conflict pits two different concrete top-level
classes against each other. The 62 conflicts are 0.05% of the 124,930 output
records.

### 2.8 Reproducibility and reconciliation with the previously released library

**Run-to-run determinism.** The delivered CLI was run twice end-to-end on the
full RepBase 31.09 input. Both runs produced

```
257cc6c118e32f1644974d7cef4ae20456dddde9ea3354f30132626d63c73299
```

and `cmp` reported the library and all five report files byte-identical, so the
conversion is deterministic.

**Difference from the library shipped with the first release.** That library
(`72f41f7532a744b82902b11e898c309562fa14c25fdb6ee52d333da3a17611ad`) was
produced before the input-validation and classification fixes documented above.
The two libraries agree on every sequence:

| Check | Result |
| --- | --- |
| Record count | 124,930 = 124,930 |
| Headers present in both | 124,872 |
| — of those, records whose sequence differs | **0** |
| Headers present in only one | 58 (class label only) |

All 58 differences are class labels, and all follow from the two rules in §2.7:

| Cause | Records |
| --- | --- |
| Identical-sequence class conflict resolved by the §2.7 rule (the earlier library used the class of the first record read) | 47 |
| Sequence classified from its own record's `KW` instead of inheriting another record's class | 11 |
| **Records whose sequence changed** | **0** |

The complete list, grouped by change:

| Change | n | Identifiers |
| --- | --- | --- |
| `SINE` → `SINE/tRNA` | 24 | `BoSB10A`, `BoSB11`, `BoSB12`, `BoSB13`, `BoSB14A`, `BoSB14B`, `BoSB15`, `BoSB1A`, `BoSB2`, `BoSB3`, `BoSB5A`, `BoSB5B`, `BoSB5C`, `BoSB5D`, `BoSB6A`, `BoSB6B`, `BoSB6C`, `BoSB6D`, `BoSB7A`, `BoSB8A`, `BoSB8B`, `BoSB9A`, `BoSB9B`, `SINE1_MT` |
| `Unknown` → `SINE/tRNA` | 14 | `B2_Mm1a`, `B2_Mm1t`, `B2_Mm2`, `B2_Rat1`, `B2_Rat2`, `B2_Rat3`, `B2_Rat4`, `B2_Rn`, `B2_Rn1`, `B2_Rn2`, `B2_RnY`, `BC1_Mm`, `ID_Rn1_dup2`, `ID_Rn2` |
| `DNA` → `DNA/TcMar` | 4 | `PIGMET_dup2`, `SONATA1_dup2`, `SONATA2_dup2`, `SONATA3` |
| `SINE/7SL` → `SINE` | 3 | `AluMacYa3_dup2`, `AluMacYb2_dup2`, `AluMacYb4_dup2` |
| `Unknown` → `SINE/7SL` | 3 | `B1_Mm`, `B1_Mus1`, `B1_Rn` |
| `LTR` → `LTR/Gypsy` | 2 | `GmGYPSY11_I`, `RAM4_I` |
| `DNA` → `DNA/MULE` | 2 | `HT1_MT`, `SHELI_MT` |
| `LTR/Copia` → `LTR/Gypsy` | 2 | `RAM13_I_MT_dup2`, `RAM13_LTR_MT_dup2` |
| `Unknown` → `DNA/hAT` | 1 | `TE-7-1_VV` |
| `SINE` → `LINE/RTE` | 1 | `TS2` |
| `SINE/tRNA` → `Unknown` | 1 | `BC1_Cp_dup2` |
| `Structural_RNA` → `Unknown` | 1 | `UHG_dup2` |

Net effect on the record-level distribution: `Unknown` 863 → 847 (−16),
`SINE` 2,915 → 2,930, `DNA` 27,739 → 27,740, `LINE` 13,511 → 13,512,
`Structural_RNA` 140 → 139; every other class is unchanged. In each case the
outcome follows the rules in §2.7 — no rule was tuned to force a hash match, and
the two libraries differ only in these 58 headers.

Runtime: 104 s wall (~100% CPU, ≈880 MB peak RSS) for the full 126,885-record
input.

---

## 3. Classification validation

### 3.1 Cross-validation against the 2018 RepeatMasker Edition

Command:

```bash
python scripts/cross_validate.py \
    --old-lib /path/to/repbase_2018.lib \
    --new-lib repbase_31.09.rm.fa
```

| Metric | Result |
| --- | --- |
| Shared names | 46,105 |
| Top-level class agreement | 45,230 / 46,105 = **98.10%** |
| Full class/subclass agreement | **60.29%** |

> **Note:** 98.10% is an *agreement rate with the 2018 RepeatMasker Edition on
> shared names*, not a classification "accuracy". There is no independent
> ground truth here; the two libraries differ partly because they make
> deliberately different modelling choices, not only because one is right.

Differences fall into four categories:

1. **Subclass granularity (majority of disagreements).** Terrier stops at the
   superfamily level, so `DNA/TcMar-Tc1` → `DNA/TcMar`, `LTR/ERV1` → `LTR/ERV`,
   `DNA/hAT-Ac` → `DNA/hAT`. This is an intentional coarsening, not an error.
2. **Genuine modern reclassification.** `LINE/Penelope` → `PLE` (307 records).
   RepeatMasker 4.x treats Penelope as its own class, so the converted library is
   **more current** than the 2018 edition here.
3. **`new = Unknown`** (222 records, 0.48%) — mostly structural RNAs
   (`4.5SRNA`, `5S`, `7SK`, `5SrRNA`) absent from the Terrier mapping.
4. **`tRNA` → `Structural_RNA`** (72 records) — the 2018 edition uses `tRNA` as a
   class name; Terrier groups it under `Structural_RNA`. Both are defensible.

### 3.2 Mapping behaviour (unit-tested)

The synthetic tests in `tests/test_conversion.py` pin down the classification
contract: case-insensitive matching, first-match-wins, `Helitron` → `RC/Helitron`
(not `DNA`), `Penelope` → `PLE`, and retention of `Unknown`.

### 3.3 Why Terrier's mapping and not `embl2rm.pl`

The alternative Cornell `embl2rm.pl` route was evaluated and rejected: it
hard-codes only four class codes and would misclassify `Helitron` as
`DNA/Helitron` (it should be `RC/Helitron`) and cannot express SINE subtypes.
See the project history for the comparison.

---

## 4. Sequence integrity

The delivered library was reconciled against the official RepBase 31.09 FASTA
release (`RepBase31.09.fasta.tar.gz`, 39 main files) by comparing the **set of
SHA-256 hashes of the distinct normalized sequences**. This is a *set*
comparison, not a multiset comparison: the two sides deliberately differ in
record count after de-duplication.

Normalization applied for the comparison (identical to the converter): uppercase,
`x→N`, `u→T`, `o→N`.

| Item | Official FASTA | This conversion |
| --- | --- | --- |
| Records | 126,885 | 124,930 |
| Distinct normalized sequences | 124,924 | **124,924** |
| Present here but not in official | — | **0** |
| Present in official but not here | 1,955 duplicate record instances (each already represented here by an identical sequence) | — |

**Conclusion:** after sequence normalization, the converted library and the
official FASTA contain exactly the same set of distinct sequences
(124,924 = 124,924). Relative to the original records, the conversion
consolidates 1,955 redundant records under the stated de-duplication rules and
loses no unique sequence; no retained sequence was altered (a real substitution,
deletion or insertion would have broken the hash match). The 1,955 official-only
records are exactly the "same ID + identical sequence" duplicates dropped by
rule (matching §2.1); in other words 1,955 is a difference in **record
instances**, not in distinct sequences.

The 124,930 output records correspond to 124,924 distinct sequences: 6 records
share a sequence with another record under a different identifier:

```
ID4_Rn#Unknown    = ID4_#SINE/tRNA
ID_Rn1#Unknown    = ID_Rn2#Unknown
B4#SINE/tRNA      = B4_ROD#Unknown
B1F1#SINE/7SL     = B1F1_mur#Unknown
B4A#SINE/tRNA     = B4A_ROD#Unknown
B1F2#SINE/7SL     = B1F2_mur#Unknown
```

These are distinct identifiers whose consensus sequences happen to be identical
(typically the same element registered once as a specific SINE and once as a
rodent entry). The identifiers do **not** conflict — de-duplication is
identifier-scoped and the output contains no duplicated identifier (§2.6) — so no
further de-duplication is applied; global sequence uniqueness is not a
requirement.

Additional structural checks:

- `makeblastdb -parse_seqids` succeeds on the full library.
- All identifiers are unique and contain only RepeatMasker-legal characters.
- Every sequence contains only IUPAC nucleotide characters.

---

## 5. RepeatMasker compatibility

### 5.1 Database build

```bash
makeblastdb -in repbase_31.09.rm.fa -dbtype nucl -parse_seqids
```

Result: **Passed** (124,930 sequences; no duplicate-`seq_id` error, confirming
de-duplication).

### 5.2 Functional test with positive controls

```bash
RepeatMasker -e rmblast -lib repbase_31.09.rm.fa -xsmall -gff \
    -dir rm_test_out test_genome.fa
```

| Control | Expected | Observed |
| --- | --- | --- |
| Gypsy-30_AnMe-I | `LTR/Gypsy` | `LTR/Gypsy` ✓ |
| L1-62_AT | `LINE/L1` | `LINE/L1` ✓ |
| Helitron-N7_AT | `RC/Helitron` | `RC/Helitron` ✓ |
| L1-3_CR (Penelope) | `PLE` | `PLE` + internal `(GT)n` Simple_repeat ✓ |
| IKIRARA1 (unclassified) | `Unknown` | `Unknown` ✓ |
| Random spacer | no hit | no hit ✓ |

The `.tbl` reported **bases masked 92.51%** — this is the expected result for a
**synthetic, positively enriched** control set (every sequence was drawn from the
library itself). **It is not a genome TE content and not a sensitivity figure.**

> The `test_genome.fa` positive-control set contains RepBase-derived sequences and
> is therefore **not** distributed with this repository.

---

## 6. Real-genome evaluation (NTN beetle, chr2 first 10 Mb)

A relative comparison was run on a real insect genome to characterise behaviour,
comparing three libraries under identical conditions. It addresses three
questions: (1) how much additional annotation RepBase 31.09 provides over the
2018 edition, in total coverage and in newly covered (non-overlapping) bases;
(2) whether that added coverage falls in classified TE categories rather than in
low-complexity or `Unknown`; and (3) the runtime cost of the larger library.

**Setup:** NTN beetle genome (409.8 Mb, 8 chromosomes), chr2 first 10 Mb (0% N);
same machine and environment (RepeatMasker 4.2.4 / RMBlast 2.17.1); identical
parameters `-xsmall -gff -pa 32`; only the library changes.

The two RepBase runs supply a custom library explicitly via `-lib`, whereas the
Dfam 4.0 run uses RepeatMasker's managed library selection (`-species "Insecta"`).
The comparison therefore characterises the effect of each **workflow** (custom
library vs. installed species library) rather than the isolated effect of
library size.

| Run | Library | Library size | Time (10 Mb) |
| --- | --- | --- | --- |
| A | RepBase 2018 RM Edition | 49,067 entries | 108 s |
| B | RepBase 31.09 (this conversion) | 124,930 entries | 252 s |
| C | Dfam 4.0 `-species "Insecta"` | 6,831 families | 54 s |

`.tbl` totals:

| Category | 2018 | 31.09 | Δ | Dfam 4.0 |
| --- | --- | --- | --- | --- |
| **Total masked bp** | 14.66% | **18.08%** | **+3.42pp** | 15.40% |
| Retroelements | 7.27% | 9.77% | +2.50pp | 5.66% |
| — LTR | 3.11% | 5.04% | +1.93pp | 2.74% |
| — LINE | 4.14% | 4.47% | +0.33pp | 2.90% |
| — Penelope (PLE) | 0.00% | 0.23% | +0.23pp | 0.01% |
| — SINE | 0.01% | 0.03% | +0.02pp | 0.01% |
| DNA transposons | 3.32% | 4.17% | +0.85pp | 2.69% |
| Rolling-circles (Helitron) | 0.28% | 0.36% | +0.08pp | 0.05% |
| Unclassified | 0.01% | 0.01% | 0.00pp | **3.30%** |
| Small RNA | 2.78% | 2.80% | +0.02pp | 2.72% |
| Simple repeats | 0.88% | 0.86% | −0.02pp | 0.88% |
| Low complexity | 0.10% | 0.10% | 0.00pp | 0.10% |

Run B was repeated with the fixed 31.09 library (the one reconciled in §2.8,
`257cc6c1…`). The resulting `.tbl` is byte-identical to the one above (SHA-256
`2d5907aa292600fdb9de9477e1aa6da2b24a4ce8f38931a3820efdc8201bbd28`; run time
261 s vs 252 s on the same machine). None of the 58 class-label changes
documented in §2.8 therefore affects any per-category percentage reported here.

Base-level overlap (31.09 vs 2018): `2018-only = 0.04%`, `31.09-only = 3.46%`,
three-way core = 1,071,066 bp.

Observations:

1. **31.09 retains the great majority of the 2018 detections.** Only 0.04% of
   the window is masked by 2018 but not by 31.09, so switching libraries keeps
   nearly all 2018 detections while adding new ones. Because that 0.04% of 2018
   detections is *not* recovered, 31.09 is not a strict superset of 2018.
2. **The added coverage is mainly attributable to classified TE categories; no
   notable low-complexity enrichment was observed.** The +3.42pp gain lands in
   LTR (+1.93pp), DNA (+0.85pp), LINE (+0.33pp), PLE (+0.23pp) and RC (+0.08pp),
   while low-complexity is 0.00pp and simple repeats −0.02pp. This shows where the
   extra coverage is *classified*; it does not by itself prove those
   classifications are biologically correct.
3. **Structural RNA is relabelled, not newly detected:** `SSU-rRNA_Hsa` and
   `LSU-rRNA_*` have identical hit counts in both libraries (e.g. 107 = 107,
   130 = 130); only the label changes (`rRNA`/`snRNA`/`tRNA` → `Structural_RNA`).
   Small RNA totals are essentially unchanged (2.78% → 2.80%). Non-TE RNA records
   are carried into masking by both libraries equally — this is not a behaviour
   introduced by 31.09.
4. **Dfam 4.0 is complementary, not additive:** 3.30pp of its 15.40% is
   `Unclassified` (curated consensus without TE class labels). Dfam also provides
   unique annotated regions that RepBase does not cover (361,797 bp masked by
   Dfam but by neither 31.09 nor 2018); this comparison does not establish the
   family origin of those regions. Results from the two sources **must not be
   added naively**; de-duplication and an explicit priority rule are required at
   the integration layer.
5. **Runtime:** 252 s per 10 Mb (2.3× the 2018 library), extrapolating to ≈2.8 h
   for a 400 Mb genome at 32 threads.

This is a **relative behavioural comparison on one 10 Mb window**: it supports
functional compatibility and annotation-behaviour assessment on real sequence,
but it does not by itself establish genome-wide annotation accuracy.

---

## 7. Known limitations

- The output is a **custom `-lib` library**, not an officially curated RepBase
  RepeatMasker Edition. It is not registered as "validated" anywhere.
- Classification is a **keyword mapping**, not sequence-based prediction.
  Records whose `KW` carries no mapped superfamily stay `Unknown`.
- **847 output records (0.68%) remain `Unknown`**; ≈180 further entries name a
  family that the mapping does not cover and could be recovered with a recorded
  synonym supplement, but this was deliberately not applied.
- **Subclass granularity is coarser** than the 2018 edition for Tc1/Mariner,
  ERV1/ERVK, hAT-Ac etc. Only the superfamily level is guaranteed.
- The library is **not restricted to any taxon**; it contains records from all
  RepBase species.
- The library **includes non-TE categories** (`Satellite`, `Structural_RNA`,
  `Unknown`, `Other`). Downstream users (e.g. gene prediction) may wish to
  annotate or exclude these separately from TE classes.
- Only the 39 main `.ref` files are converted; the `appendix/` records are
  excluded by design.
- Genome-wide sensitivity/specificity is **not** established by this report.
- The runtime penalty scales with library size (≈2.3× relative to the 2018
  edition in the tested window).

---

## 8. Reproduction

1. Obtain RepBase 31.09 under the applicable GIRI license and extract the 39
   main `.ref` files (keep the `appendix/` subdirectory out of the input
   directory).
2. Run the converter (§2) and compare the printed summary to §2.1–§2.3.
3. Cross-check integrity against the official RepBase FASTA release by
   normalizing sequences (uppercase, `x→N`, `u→T`, `o→N`) and comparing the
   set of SHA-256 hashes of the distinct sequences; expect 124,924 distinct
   sequences on both sides.
4. Optionally cross-validate classification against a reference library with
   `scripts/cross_validate.py` (§3.1).
5. Build the BLAST database (§5.1) and run the RepeatMasker functional test
   (§5.2) with your own positive controls.
6. Run the synthetic unit tests: `python tests/test_conversion.py`.

No RepBase data is distributed with this repository, so steps 1 and 5 require
licensed input obtained by the user.
