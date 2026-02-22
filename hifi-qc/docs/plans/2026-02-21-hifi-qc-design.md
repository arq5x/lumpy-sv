# hifi-qc Design Document

**Date:** 2026-02-21
**Status:** Approved

## Overview

hifi-qc is a command-line tool for quality control of PacBio HiFi DNA sequencing data. It takes one or more BAM/CRAM files as input and generates a single self-contained interactive HTML report in the style of MultiQC, with all samples overlaid on each chart for easy comparison.

## Tech Stack

- **Rust core** compiled as a Python extension module via PyO3/maturin
- **noodles** for BAM/CRAM I/O (pure Rust, no C dependencies)
- **rayon** for parallelization of coverage computation
- **Python CLI** using click for the user-facing interface
- **Jinja2 + Plotly.js** for interactive HTML report generation
- Distributed as a single `pip install hifi-qc` package

## Architecture

```
┌─────────────────────────────────────┐
│  Python CLI (click)                 │  User-facing CLI
│  hifi-qc *.bam -o report.html      │
├─────────────────────────────────────┤
│  Python Report Engine               │  Jinja2 + Plotly.js
│  Renders HTML from structured data  │
├─────────────────────────────────────┤
│  Rust Core (PyO3 extension)         │  Performance-critical
│  - Single-pass BAM/CRAM processing  │
│  - D4-style coverage engine         │
│  - Per-read metric extraction       │
│  Returns: Python dicts/lists        │
└─────────────────────────────────────┘
```

## QC Modules (11 report sections)

| # | Module | Rust Output | Plot Type |
|---|--------|-------------|-----------|
| 1 | Summary Table | Basic stats (reads, bases, N50, mean length, mean Q, mean coverage) | HTML table |
| 2 | Read Length Distribution | Length histogram bins + N50/mean/median | Overlaid histograms |
| 3 | Q-Score Distribution | Per-read Q-score histogram bins | Overlaid histograms |
| 4 | Read Length vs. Quality | (length, quality) pairs, subsampled for plotting | Scatter/heatmap |
| 5 | HiFi Passes | np tag distribution + rq tag distribution | Dual histogram |
| 6 | Genome-wide Coverage | Mean, median, stdev + coverage distribution | Histogram |
| 7 | Per-Chromosome Coverage | Mean coverage per reference sequence | Grouped bar chart |
| 8 | Coverage Uniformity | Genome fraction at >= X coverage curve | Line plot |
| 9 | GC Bias | Observed GC content vs. expected + coverage vs. GC | Dual line plot |
| 10 | Mismatch Analysis | Substitution type counts (A>G, C>T, etc.) | Stacked bar chart |
| 11 | Indel Analysis | Insertion/deletion size distribution, separate rates | Histogram + summary |

## Single-Pass Processing

All metrics are collected in a single sequential pass through each BAM/CRAM file:

```
For each BAM/CRAM:
  1. Read header -> allocate depth arrays per chromosome
  2. Iterate all records (single pass):
     - Record length, Q-score, MAPQ, np/rq tags -> accumulators
     - Walk CIGAR -> increment depth array + count mismatches/indels
     - Compute read GC content
  3. After iteration: summarize depth arrays -> coverage stats
```

### Coverage Algorithm (D4-inspired)

The coverage computation uses principles from the D4 format:

- **Pre-allocated depth arrays**: One `Vec<u32>` per chromosome, sized from BAM header. Memory proportional to genome length, not read depth.
- **CIGAR-aware depth increment**: Walking each read's CIGAR ops, incrementing depth only for M/=/X operations. Properly handles soft clips (excluded), deletions (zero depth), and insertions (no reference consumption).
- **Single pass**: Depth arrays are updated as reads are encountered during the sequential iteration. No separate pass needed.
- **Post-pass summarization**: After all reads processed, depth arrays are summarized in parallel (rayon) for per-chromosome means, genome-wide stats, coverage histogram, and GC-stratified coverage.

### GC Bias Computation

Requires a reference FASTA (optional `--reference` flag):
- During the single pass, each read's GC content is computed from its sequence.
- After the pass, reference GC content is computed in windows (e.g., 100bp) from the FASTA.
- Coverage is stratified by GC bin to produce coverage-vs-GC curves.
- If no reference is provided, GC bias section is skipped.

## Sampling Mode

For quick QC on large files, a fraction-based sampling mode is available:

```bash
hifi-qc *.bam --sample-fraction 0.10 -o report.html
```

- Each read has a `sample_fraction` probability of being processed (probabilistic, streaming).
- Uses a seedable RNG (`--seed` flag for reproducibility).
- **Coverage sections are skipped** when sampling is active (fractional depth is not meaningful).
- All other sections (length, quality, passes, errors, GC content, MAPQ) work from sampled reads.
- Report header notes the sampling fraction used.

## CLI Interface

```bash
# Basic usage
hifi-qc sample1.bam sample2.bam -o report.html

# With options
hifi-qc *.bam \
  --output report.html \
  --threads 8 \
  --reference ref.fa \
  --sample-names s1,s2 \
  --sample-fraction 0.10 \
  --seed 42
```

### Arguments

| Flag | Default | Description |
|------|---------|-------------|
| `BAM/CRAM` (positional) | required | One or more input alignment files |
| `--output, -o` | `hifi_qc_report.html` | Output HTML report path |
| `--threads, -t` | all available | Number of threads for parallel operations |
| `--reference, -r` | none | Reference FASTA (needed for GC bias) |
| `--sample-names` | from BAM RG | Comma-separated sample names |
| `--sample-fraction` | 1.0 (all reads) | Fraction of reads to sample (0.0-1.0) |
| `--seed` | none | RNG seed for reproducible sampling |

## Multi-Sample Handling

- All samples overlaid on each chart with a consistent color palette.
- Summary table has one row per sample.
- Legend identifies each sample by name (from `--sample-names` or BAM read group).
- Works with any reference genome (chromosome names/lengths from BAM header).

## Report Output

- Single self-contained HTML file (all JS/CSS inlined).
- Plotly.js embedded for interactive charts (zoom, hover, pan).
- Navigation sidebar with links to each section.
- General Statistics table at top.
- Responsive layout.

## Project Structure

```
hifi-qc/
├── Cargo.toml
├── pyproject.toml
├── src/
│   ├── lib.rs                 # PyO3 module entry point
│   ├── bam_reader.rs          # Shared BAM/CRAM reading (noodles)
│   ├── coverage.rs            # Depth array + coverage stats
│   ├── fragments.rs           # Fragment length extraction
│   ├── passes.rs              # HiFi pass/quality tag extraction
│   ├── errors.rs              # Mismatch + indel analysis (CIGAR/MD)
│   ├── gc_bias.rs             # GC content computation
│   ├── mapq.rs                # Mapping quality extraction
│   └── qscore.rs              # Q-score extraction
├── python/
│   └── hifi_qc/
│       ├── __init__.py
│       ├── cli.py             # Click CLI
│       ├── report.py          # Report generator
│       └── templates/
│           └── report.html.j2 # Jinja2 template
├── tests/
│   ├── rust/
│   └── python/
└── docs/
    └── plans/
```

## Dependencies

### Rust (Cargo.toml)
- noodles (BAM/CRAM/FASTA reading)
- pyo3 (Python bindings)
- rayon (parallelism)
- rand (sampling RNG)

### Python (pyproject.toml)
- maturin (build system)
- click (CLI)
- jinja2 (HTML templating)
- plotly (chart data serialization, optional)
