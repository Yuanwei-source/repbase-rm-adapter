# RepBase-to-RepeatMasker Adapter

**English** · [简体中文](README.zh-CN.md)

Convert RepBase EMBL records into a RepeatMasker-compatible FASTA library with
reproducible classification mapping, sequence normalization, duplicate handling,
and traceable conversion reports.

**Tested with:** RepBase 31.09 · RepeatMasker 4.2.4 · RMBlast 2.17.1

> This is an independent conversion utility, not an official RepBase RepeatMasker
> Edition release. RepBase sequence data are **not** distributed with this
> repository.

## Overview

Recent RepBase releases provide repetitive-element sequences and annotations in
EMBL format. RepeatMasker can align against a custom FASTA library directly via
`-lib`, but its per-class annotation and summary output depend on the
`>id#class/subclass` FASTA header convention; a raw RepBase FASTA therefore
produces hits that RepeatMasker cannot classify or summarise correctly.

This utility provides:

- Conversion of RepBase EMBL `.ref` records into `>sequence_id#class/subclass`
  FASTA format.
- Classification using a RepBase-to-RepeatMasker mapping adapted from the
  [Terrier](https://github.com/rbturnbull/terrier) project.
- Case-insensitive keyword matching, retaining `Unknown` when classification is
  unresolved.
- Duplicate identifier resolution with traceable renaming.
- Normalization of non-IUPAC sequence characters.
- Generation of local conversion and audit reports.

The original RepBase files are never modified.

## Requirements

- **Python ≥ 3.11** (uses the standard-library `tomllib`; tested with 3.13.13).
  No third-party Python packages are required.
- RepBase EMBL records obtained under the appropriate GIRI license.
- [RepeatMasker](https://www.repeatmasker.org/) (tested with 4.2.4).
- [RMBlast](https://www.repeatmasker.org/RMBlast.html) for repeat annotation
  (tested with 2.17.1).
- BLAST+ `makeblastdb` for optional database-index validation.

No Terrier neural-network inference or REPET installation is required. The
converter uses the **classification mapping table**, not the Terrier prediction
model.

## Usage

### 1. Obtain RepBase EMBL data

Download and extract the required RepBase release from the official distribution
source. Keep the extracted source files in a directory outside this repository.

The conversion tested here uses the 39 main `.ref` files from RepBase 31.09,
**excluding** the `appendix/` subdirectory (which shares many record IDs with the
main files). The converter reads only top-level `*.ref` files.

### 2. Convert the database

```bash
python repbase_to_rm.py \
    --input-dir /path/to/RepBase31.09.embl \
    --mapping   mappings/repbase-to-repeatmasker.toml \
    --output    output/repbase_31.09.rm.fa \
    --out-dir   output/reports \
    --allow-length-mismatch
```

`--allow-length-mismatch` is required for the RepBase 31.09 release: five records
in `humsub.ref` (`AluY`, `AluYa5`, `AluYa8`, `AluYb8`, `AluYb9`) declare an `SQ`
base-pair count that disagrees with the length of their sequence lines. Without
the flag the converter aborts on such a record instead of silently writing a
truncated or padded sequence; with it, every mismatch is listed individually in
the summary. The flag exists to accommodate known anomalies in a specific source
release; it is **not** a recommended default. For any other RepBase version, run
without it first and inspect each reported mismatch before deciding to accept it.

Audit reports (`id_mapping.tsv`, `unmapped_records.tsv`,
`base_normalization.tsv`, `duplicate_report.tsv`; plus `class_conflicts.tsv` when
records sharing an identical sequence carry different classifications) are
written next to the output library, or to `--out-dir` when given. Run
`python repbase_to_rm.py --help` for all options.

### 3. Validate the FASTA library

```bash
makeblastdb -in output/repbase_31.09.rm.fa -dbtype nucl -parse_seqids
```

The completed database should contain unique sequence identifiers and valid
nucleotide sequences.

Optionally cross-validate the classification against a reference library (for
example the official 2018 RepeatMasker Edition):

```bash
python scripts/cross_validate.py \
    --old-lib /path/to/repbase_2018.lib \
    --new-lib output/repbase_31.09.rm.fa
```

### 4. Run RepeatMasker

```bash
RepeatMasker -e rmblast \
    -lib output/repbase_31.09.rm.fa \
    -xsmall -gff \
    -dir rm_results \
    genome.fa
```

The converted library is supplied explicitly through `-lib`; it does not replace
RepeatMasker's installed reference libraries. Users may filter the resulting
annotations according to their organism, research question and repeat
categories.

### 5. Tests

Synthetic-data tests (no RepBase sequences required, standard library only):

```bash
python tests/test_conversion.py
# or: python -m unittest discover -s tests -v
```

## Conversion strategy

**Classification.** The first case-insensitive matching keyword in each record's
`KW` field is mapped to the RepeatMasker classification hierarchy. Entries
without a recognized classification remain `Unknown`; they are not assigned a
category by guessing.

**Classification conflicts.** When one identifier covers several records whose
sequences are identical after normalization but whose keyword fields map to
different classes, the converter resolves the conflict deterministically — a
concrete class beats `Unknown`, a class with a subclass beats a bare top-level
class, and otherwise the first record in input order wins — and records every
decision in `class_conflicts.tsv`. Conflicts are never resolved silently.

**Duplicate identifiers.** Identical normalized sequences sharing an identifier
are consolidated. Distinct sequences sharing an identifier receive unique
suffixes (`_dup2`, `_dup3`, …). Duplicate comparisons are performed after
sequence normalization.

**Sequence normalization.** Non-IUPAC `X` and `O` characters are converted to
`N`, and `U` is converted to `T`. Valid IUPAC ambiguity characters are preserved.

**Identifier sanitization.** Characters that are illegal in RepeatMasker FASTA
identifiers (`:`, `(`, `)`, `@`) are replaced with `_`.

**Traceability.** The conversion retains local reports linking original
identifiers, output identifiers, classification decisions and
sequence-normalization events.

## Validation

The adapter was evaluated using RepBase 31.09 and RepeatMasker 4.2.4 on
2026-10-09.

| Metric | Result |
| --- | --- |
| Input EMBL records | 126,885 |
| Output FASTA records | 124,930 |
| Classified records | 124,083 |
| Unclassified records | 847 (0.68%) |
| Classification coverage | 99.32% |
| Top-level class agreement with the 2018 RepeatMasker Edition (shared names) | 98.10% |
| `makeblastdb -parse_seqids` | Passed |
| RepeatMasker functional test | Passed |
| Distinct normalized sequences vs official FASTA | 124,924 = 124,924 (identical set) |
| Real-genome comparison (NTN beetle, chr2 first 10 Mb) | 31.09 = 18.08% masked vs 2018 = 14.66%; +3.42pp mainly in classified TE |

Functional testing confirmed annotation and classification output for Gypsy, L1,
Helitron, Penelope and an unclassified repeat. A random spacer control produced
no reported repeat match.

These tests establish conversion compatibility and basic annotation
functionality. They do **not** establish genome-wide annotation sensitivity or
specificity, nor superiority over the curated 2018 RepeatMasker Edition.

Full validation details, including the real-genome comparison, are provided in
[`docs/VALIDATION.md`](docs/VALIDATION.md).

## Limitations

- The output is a custom library, not an officially curated RepBase RepeatMasker
  Edition.
- Classification is based on an explicit keyword mapping, not sequence-based
  prediction.
- Some subtype distinctions in the 2018 RepeatMasker Edition are represented at
  broader superfamily levels (e.g. `DNA/TcMar-Tc1` → `DNA/TcMar`).
- The full conversion library is not automatically restricted to a target
  taxonomic group.
- The converted library includes non-TE repeat categories (satellites,
  structural RNA); users should select appropriate categories for downstream
  applications.
- 847 output records (0.68%) remain `Unknown` because their `KW` field carries no
  mapped superfamily; they are retained, not discarded.

## Data availability and licensing

RepBase is maintained by the Genetic Information Research Institute (GIRI).
RepBase has historically been distributed under a restrictive data-use agreement,
in which the academic-user agreement restricts redistribution of the database and
of materials derived from it outside the research group.

2026 open-licence development. In June 2026, GIRI and the Dfam team announced
that RepBase is being integrated into a unified open database framework, and
decided to release the complete RepBase dataset under CC0 (Creative Commons
Zero) (Kojima et al., 2026, Mobile DNA, 17:17;
[DOI: 10.1186/s13100-026-00409-9](https://doi.org/10.1186/s13100-026-00409-9)).

The integration will be carried out in stages. The announcement of the CC0 plan
should not be equated with every historical or current RepBase release having
already completed its licence change. When using or redistributing the data of a
particular version, the licence terms explicitly provided by the corresponding
data provider for that version shall govern.

This repository does not distribute RepBase data, including the raw EMBL
records, the converted consensus sequences, and test data derived from RepBase
sequences. Users can obtain the required version from the
[official GIRI download page](https://www.girinst.org/server/RepBase/) and must
comply with the corresponding licence terms.

The classification mapping is adapted from the Terrier project (Apache-2.0); the
relevant sources and licences are given in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

The original code in this repository is released under the MIT License; see
[LICENSE](LICENSE). That licence applies only to this project's original code
and does not alter the respective licence conditions of RepBase data or
third-party materials.

## References and acknowledgments

- **RepBase:** Genetic Information Research Institute — <https://www.girinst.org/>
  (downloads: <https://www.girinst.org/server/RepBase/>)
- **RepeatMasker:** <https://www.repeatmasker.org/>
- **Terrier:** <https://github.com/rbturnbull/terrier>
- **Terrier publication:** Turnbull et al. (2025). *Terrier: a deep learning
  repeat classifier*. Briefings in Bioinformatics, 26(4), bbaf442.
  <https://doi.org/10.1093/bib/bbaf442>
- **RepBase–Dfam unification (2026):** Kojima, K. K., Smit, A. F. A., Bao, W.,
  Kohany, O., Kojima, N. F., Jurka, T., Hubley, R., & Wheeler, T. J. *Unifying
  Repbase and Dfam: a new open foundation for transposable element research*.
  Mobile DNA, 17, Article 17 (2026). <https://doi.org/10.1186/s13100-026-00409-9>
  (PMID 42343465). GIRI statement:
  <http://www.girinst.org/repbase/repbase_and_dfam.html>

## Project status

The conversion workflow has passed format validation and RepeatMasker functional
testing with RepBase 31.09. It is distributed as an independently developed,
reproducible conversion utility. Downstream genome-annotation performance should
be assessed separately for each application.

**Release v1.1.0 — data-integrity and classification improvements, validated on
the full RepBase 31.09 library; the validated conversion and functional
compatibility with RepeatMasker 4.2.4 established in v1.0.0 are unchanged.**
"Validated" here means validated by this project's own tests (conversion
integrity, RepeatMasker compatibility, and a real-genome behavioural
comparison); it does **not** mean GIRI certification or a genome-wide
annotation-accuracy guarantee.
