# hifi-qc

Quality control for PacBio HiFi DNA sequencing data. Generates interactive HTML reports in the style of [MultiQC](https://multiqc.info/).

Takes one or more BAM/CRAM files as input and produces a single self-contained HTML report with all samples overlaid on each chart for easy comparison.

## Features

- **Single-pass processing** -- all metrics collected in one pass through each BAM/CRAM
- **D4-inspired coverage engine** -- pre-allocated depth arrays, CIGAR-aware, parallel summarization via rayon
- **11 QC modules** -- read length, Q-score, HiFi passes, coverage, GC bias, mapping quality, mismatches, indels
- **Interactive charts** -- Plotly.js with zoom, hover, and pan
- **Multi-sample support** -- overlaid plots with color-coded samples
- **Sampling mode** -- fast QC on large files by processing a fraction of reads
- **BAM and CRAM support** -- via noodles (pure Rust, no htslib dependency)

## Report Sections

| Section | Description |
|---------|-------------|
| General Statistics | Summary table: reads, bases, N50, mean Q, coverage, mapped % |
| Read Length Distribution | Histogram of HiFi read lengths |
| Q-Score Distribution | Per-read mean quality score histogram |
| Read Length vs Quality | Scatter plot revealing length-quality correlation |
| HiFi Passes | Distribution of CCS pass counts (np tag) |
| Genome-wide Coverage | Per-base depth distribution |
| Per-Chromosome Coverage | Mean coverage per reference sequence |
| Coverage Uniformity | Genome fraction at >= X coverage |
| GC Bias | Coverage vs GC content (requires `--reference`) |
| Mapping Quality | MAPQ score distribution |
| Mismatches & Indels | Insertion/deletion size distributions |

## Installation

### Prerequisites

- Python >= 3.10
- Rust toolchain (install via [rustup](https://rustup.rs/))

### From source

```bash
git clone https://github.com/arq5x/hifi-qc.git
cd hifi-qc
python -m venv .venv
source .venv/bin/activate
pip install maturin
maturin develop --release
```

This compiles the Rust core and installs `hifi-qc` as a command-line tool in your virtual environment.

## Usage

```
Usage: hifi-qc [OPTIONS] BAM_FILES...

  Quality control for PacBio HiFi sequencing data.

Options:
  -o, --output TEXT        Output HTML report path
  -r, --reference PATH     Reference FASTA for GC bias
  -t, --threads INTEGER    Number of threads
  --sample-fraction FLOAT  Fraction of reads to sample (0.0-1.0)
  --seed INTEGER           RNG seed for reproducible sampling
  --sample-names TEXT      Comma-separated sample names
  --help                   Show this message and exit.
```

### Basic usage

```bash
hifi-qc sample.bam -o report.html
```

### Multiple samples

All samples are overlaid on each chart with different colors:

```bash
hifi-qc sample1.bam sample2.bam sample3.bam -o comparison.html
```

### With reference (enables GC bias analysis)

Provide an indexed FASTA (`.fai` must exist alongside the `.fa`):

```bash
hifi-qc sample.bam -o report.html -r reference.fa
```

### Custom sample names

Override the names shown in the report (default: filename stem or BAM read group):

```bash
hifi-qc s1.bam s2.bam -o report.html --sample-names "Patient_A,Patient_B"
```

### Sampling mode (fast QC)

Process only a fraction of reads for a quick look at large files. Coverage sections are skipped in this mode since fractional depth is not meaningful:

```bash
hifi-qc large.bam -o quick_report.html --sample-fraction 0.10 --seed 42
```

### Multi-threaded

Coverage summarization is parallelized across chromosomes:

```bash
hifi-qc sample.bam -o report.html -t 8
```

### Example with test data

The repository includes a small synthetic BAM for testing:

```bash
# From the project root, with the venv activated:
hifi-qc tests/data/test.bam -o test_report.html
```

With GC bias analysis:

```bash
hifi-qc tests/data/test.bam -o test_report.html -r tests/data/test_ref.fa
```

Open `test_report.html` in a browser to see the interactive report.

## Architecture

```
┌─────────────────────────────────────┐
│  Python CLI (click)                 │  User-facing interface
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

The Rust core does a single sequential pass through each BAM/CRAM file, simultaneously collecting all per-read metrics (length, quality, MAPQ, HiFi tags, GC content) and updating depth arrays via CIGAR walking. After the pass, coverage statistics are summarized in parallel across chromosomes using rayon. Results are returned as Python dictionaries to the report engine.

### Project structure

```
hifi-qc/
├── Cargo.toml                          # Rust dependencies
├── pyproject.toml                      # Python package config (maturin)
├── src/
│   ├── lib.rs                          # PyO3 entry point
│   ├── bam_reader.rs                   # BAM/CRAM reading (noodles)
│   ├── metrics.rs                      # Per-read metric accumulator
│   ├── coverage.rs                     # Parallel coverage summarization
│   └── gc_bias.rs                      # GC bias computation
├── python/
│   └── hifi_qc/
│       ├── cli.py                      # Click CLI
│       ├── report.py                   # Report generator
│       └── templates/
│           └── report.html.j2          # HTML report template
└── tests/
    ├── create_test_bam.py              # Synthetic test BAM generator
    ├── test_e2e.py                     # End-to-end tests
    └── data/                           # Test BAM and reference
```

## Dependencies

### Rust
- [noodles](https://github.com/zaeleus/noodles) 0.88 -- BAM/CRAM/FASTA I/O (pure Rust)
- [PyO3](https://pyo3.rs/) 0.25 -- Python bindings
- [rayon](https://github.com/rayon-rs/rayon) 1.10 -- parallelism
- [rand](https://github.com/rust-random/rand) 0.9 -- sampling RNG

### Python
- [click](https://click.palletsprojects.com/) >= 8.0 -- CLI framework
- [Jinja2](https://jinja.palletsprojects.com/) >= 3.0 -- HTML templating
- [Plotly.js](https://plotly.com/javascript/) 2.35 -- interactive charts (loaded from CDN)

## Running tests

```bash
# Rust unit tests (31 tests)
cargo test --lib

# End-to-end tests (requires pysam for test data generation)
pip install pysam
python tests/create_test_bam.py   # generate test data (already included)
python tests/test_e2e.py
```

## License

TBD
