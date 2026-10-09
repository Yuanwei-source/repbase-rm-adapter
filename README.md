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
EMBL format. These records must be converted into RepeatMasker-compatible FASTA
headers before they can be used with the `-lib` option.

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
    --out-dir   output/reports
```

Audit reports (`id_mapping.tsv`, `unmapped_records.tsv`,
`base_normalization.tsv`, `duplicate_report.tsv`) are written next to the output
library, or to `--out-dir` when given. Run `python repbase_to_rm.py --help` for
all options.

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
| Classified records | 124,067 |
| Unclassified records | 863 (0.69%) |
| Classification coverage | 99.31% |
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
- 863 output records (0.69%) remain `Unknown` because their `KW` field carries no
  mapped superfamily; they are retained, not discarded.

## Data availability and licensing

RepBase is maintained by the Genetic Information Research Institute (GIRI) and
is subject to its data-use and redistribution terms. GIRI's academic-user terms
restrict redistribution of the database, its components and derived materials.

> **2026 update — RepBase moving to CC0.** In June 2026, GIRI announced that
> RepBase will be released under a **CC0 public-domain license** as part of an
> effort to unify RepBase and Dfam into a single open framework (Storer et al.,
> *Mobile DNA* 17:17, 2026; DOI
> [10.1186/s13100-026-00409-9](https://doi.org/10.1186/s13100-026-00409-9)).
> This is a **multi-release transition**: the announcement does not mean that
> Dfam 4.0 already incorporates RepBase 31.09 in full, nor that every current
> RepBase release is redistributable today. Until the transition covers the
> specific release and date you use, obtain RepBase from GIRI under the terms
> that apply to it. This repository continues to distribute no RepBase data.

This repository contains conversion software and documentation, but does **not**
distribute RepBase source records, converted consensus sequences, or
sequence-derived test data. Users must obtain RepBase independently under the
applicable license.

The classification mapping is derived from the Terrier project
(Apache-2.0). See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

Original code in this repository is released under the MIT License; see
[`LICENSE`](LICENSE). That license does not extend to RepBase data.

## References and acknowledgments

- **RepBase:** Genetic Information Research Institute — <https://www.girinst.org/>
- **RepeatMasker:** <https://www.repeatmasker.org/>
- **Terrier:** <https://github.com/rbturnbull/terrier>
- **Terrier publication:** Turnbull et al. (2025). *Terrier: a deep learning
  repeat classifier*. Briefings in Bioinformatics, 26(4), bbaf442.
  <https://doi.org/10.1093/bib/bbaf442>
- **RepBase–Dfam unification (2026):** Storer, J. M., Hubley, R. M., Rosen, J.
  B., Wheeler, T. J., & Smit, A. F. A. *Unifying Repbase and Dfam: a new open
  foundation for transposable element research*. Mobile DNA, 17, 17.
  <https://doi.org/10.1186/s13100-026-00409-9> (PMID 42343465). GIRI statement:
  <http://www.girinst.org/repbase/repbase_and_dfam.html>

## Project status

The conversion workflow has passed format validation and RepeatMasker functional
testing with RepBase 31.09. It is distributed as an independently developed,
reproducible conversion utility. Downstream genome-annotation performance should
be assessed separately for each application.

**Release v1.0.0 — validated conversion and functional compatibility with
RepBase 31.09 and RepeatMasker 4.2.4.** "Validated" here means validated by this
project's own tests (conversion integrity, RepeatMasker compatibility, and a
real-genome behavioural comparison); it does **not** mean GIRI certification or a
genome-wide annotation-accuracy guarantee.
