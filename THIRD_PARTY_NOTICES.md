# Third-Party Notices

This project includes or derives from third-party material. The original code in
this repository is licensed under the MIT License (see `LICENSE`). The notices
below cover third-party components.

## 1. Terrier classification mapping (included)

`mappings/repbase-to-repeatmasker.toml` is redistributed **unmodified** from:

- **Project:** Terrier — <https://github.com/rbturnbull/terrier>
- **File:** `terrier/data/repbase-to-repeatmasker.toml`
- **Revision:** `main` @ `dbbc35972f43ede1b863f1fade1e10545a9f020a` (2026-10-02)
- **License:** Apache License 2.0
- **SHA-256 (this copy):** `b45e04e4c1107a9ffee15debe361713fedfdce4c5ae27ad79da631ca7818e1a1`
- **Modifications:** none (byte-identical to the upstream file)

The full Apache-2.0 license text is provided at `LICENSES/Apache-2.0.txt`.
The upstream project has no `NOTICE` file. Per Apache-2.0 §4, this file retains
the original license and attribution, and states that no modifications were made.

### What this project uses

This project uses Terrier's **classification mapping table only**. It does **not**
run Terrier's neural-network prediction model, and it does not require Terrier or
its dependencies to be installed.

### Citation

Turnbull et al. (2025). *Terrier: a deep learning repeat classifier*.
Briefings in Bioinformatics, 26(4), bbaf442.
<https://doi.org/10.1093/bib/bbaf442>

## 2. RepBase (NOT included)

RepBase is maintained by the Genetic Information Research Institute (GIRI) —
<https://www.girinst.org/>. RepBase has historically been distributed under a
restrictive data-use agreement; GIRI's academic-user terms restrict
redistribution of the database and of derived materials outside the research
group.

The RepBase distributions examined for this project contain **no per-release
licence file**: the RepBase 31.09 EMBL and FASTA archives ship only `.ref`
records, and the 2018 RepeatMasker Edition source ships a format/compatibility
`README.RMRBSeqs` with no licence text. The applicable terms are therefore those
GIRI sets for the release and account under which RepBase is obtained.

**Status of the announced CC0 transition.** In June 2026, GIRI and the Dfam team
announced that RepBase is being integrated into a unified open database
framework, and decided to release the complete RepBase dataset under CC0
(Kojima et al., 2026, *Mobile DNA*, 17, Article 17;
DOI 10.1186/s13100-026-00409-9). The integration will be carried out in stages.
The announcement of the CC0 plan should **not** be equated with every historical
or current RepBase release having already completed its licence change, and this
project therefore does **not** treat any RepBase release as CC0; the
conservative position below is retained. When using or redistributing the data
of a particular version, the licence terms explicitly provided by the
corresponding data provider for that version shall govern (downloads:
<https://www.girinst.org/server/RepBase/>).

**No RepBase records, consensus sequences or sequence-derived test data are
distributed in this repository.** Users must obtain RepBase independently under
the applicable license. Conversion outputs (the `-lib` FASTA library, audit TSVs
and any RepBase-derived test data) must be kept local and are excluded by
`.gitignore`. The MIT License covering this repository's original code applies
only to that code and does not alter the licence conditions of RepBase data or
third-party materials.

## 3. RepeatMasker / RMBlast (NOT included)

[RepeatMasker](https://www.repeatmasker.org/) and
[RMBlast](https://www.repeatmasker.org/RMBlast.html) are external tools used for
validation. They are not distributed here and are governed by their own
licenses. See the RepeatMasker website for current terms.
