# hifi-qc Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a PacBio HiFi sequencing QC tool that generates MultiQC-style interactive HTML reports from BAM/CRAM files.

**Architecture:** Rust core (PyO3 extension via maturin) handles all BAM/CRAM parsing and metric computation in a single pass per file. Python layer provides the CLI (click) and renders a self-contained HTML report (Jinja2 + Plotly.js). The Rust core returns structured Python dicts/lists to the Python layer.

**Tech Stack:** Rust (noodles, pyo3, rayon, rand), Python (click, jinja2), Plotly.js (embedded in HTML)

---

## Task 1: Project scaffolding and build system

**Files:**
- Create: `Cargo.toml`
- Create: `pyproject.toml`
- Create: `src/lib.rs`
- Create: `python/hifi_qc/__init__.py`
- Create: `python/hifi_qc/cli.py`

**Step 1: Create Cargo.toml**

```toml
[package]
name = "hifi-qc"
version = "0.1.0"
edition = "2021"

[lib]
name = "hifi_qc"
crate-type = ["cdylib"]

[dependencies]
pyo3 = { version = "0.23", features = ["extension-module"] }
noodles = { version = "0.88", features = ["bam", "sam", "cram", "fasta", "core", "bgzf"] }
rayon = "1.10"
rand = "0.9"
```

**Step 2: Create pyproject.toml**

```toml
[build-system]
requires = ["maturin>=1.0,<2.0"]
build-backend = "maturin"

[project]
name = "hifi-qc"
version = "0.1.0"
description = "Quality control for PacBio HiFi sequencing data"
requires-python = ">=3.10"
dependencies = [
    "click>=8.0",
    "jinja2>=3.0",
]

[project.scripts]
hifi-qc = "hifi_qc.cli:main"

[tool.maturin]
features = ["pyo3/extension-module"]
python-source = "python"
module-name = "hifi_qc._core"
```

**Step 3: Create minimal src/lib.rs**

```rust
use pyo3::prelude::*;

/// Placeholder: process BAM/CRAM files and return QC metrics.
#[pyfunction]
#[pyo3(signature = (bam_paths, reference=None, threads=None, sample_fraction=None, seed=None))]
fn process_bam_files(
    _bam_paths: Vec<String>,
    _reference: Option<String>,
    _threads: Option<usize>,
    _sample_fraction: Option<f64>,
    _seed: Option<u64>,
) -> PyResult<Vec<PyObject>> {
    Ok(vec![])
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(process_bam_files, m)?)?;
    Ok(())
}
```

**Step 4: Create python/hifi_qc/__init__.py**

```python
"""hifi-qc: Quality control for PacBio HiFi sequencing data."""
```

**Step 5: Create minimal python/hifi_qc/cli.py**

```python
import click


@click.command()
@click.argument("bam_files", nargs=-1, required=True, type=click.Path(exists=True))
@click.option("-o", "--output", default="hifi_qc_report.html", help="Output HTML report path")
@click.option("-r", "--reference", default=None, type=click.Path(exists=True), help="Reference FASTA for GC bias")
@click.option("-t", "--threads", default=None, type=int, help="Number of threads")
@click.option("--sample-fraction", default=None, type=float, help="Fraction of reads to sample (0.0-1.0)")
@click.option("--seed", default=None, type=int, help="RNG seed for reproducible sampling")
@click.option("--sample-names", default=None, type=str, help="Comma-separated sample names")
def main(bam_files, output, reference, threads, sample_fraction, seed, sample_names):
    """Quality control for PacBio HiFi sequencing data."""
    from hifi_qc._core import process_bam_files

    results = process_bam_files(
        list(bam_files),
        reference=reference,
        threads=threads,
        sample_fraction=sample_fraction,
        seed=seed,
    )
    click.echo(f"Processed {len(bam_files)} file(s). Results: {len(results)} sample(s).")
    click.echo(f"Report would be written to: {output}")
```

**Step 6: Install maturin and verify the build compiles**

Run: `pip install maturin && cd /Users/arq5x/src/quinlan/hifi-qc && maturin develop`
Expected: Build succeeds, `hifi-qc --help` shows the CLI options.

**Step 7: Verify CLI runs**

Run: `cd /Users/arq5x/src/quinlan/hifi-qc && hifi-qc --help`
Expected: Shows help with BAM_FILES argument and all options.

**Step 8: Commit**

```bash
git add Cargo.toml pyproject.toml src/lib.rs python/hifi_qc/__init__.py python/hifi_qc/cli.py
git commit -m "feat: project scaffolding with PyO3/maturin build system"
```

---

## Task 2: BAM/CRAM reader module

**Files:**
- Create: `src/bam_reader.rs`
- Modify: `src/lib.rs`

This module provides a unified interface to iterate records from BAM or CRAM files, extracting the header (reference sequences/lengths) and yielding records.

**Step 1: Create src/bam_reader.rs**

This module defines:
- `FileHeader`: extracted reference sequence names and lengths
- `ReadRecord`: a struct holding all per-read fields we need (extracted during iteration so we don't hold borrowed BAM record data)
- `read_bam_file()` and `read_cram_file()` functions that iterate records and call a closure for each

```rust
use noodles::bam;
use noodles::cram;
use noodles::fasta;
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use std::fs::File;
use std::io;
use std::path::Path;

/// Extracted header info.
pub struct FileHeader {
    pub reference_names: Vec<String>,
    pub reference_lengths: Vec<usize>,
    pub sample_name: Option<String>,
}

/// All per-read data extracted during the single pass.
pub struct ReadRecord {
    pub sequence_length: usize,
    pub quality_scores: Vec<u8>,
    pub mapping_quality: u8,
    pub reference_id: Option<usize>,
    pub alignment_start: Option<usize>, // 0-based
    pub cigar_ops: Vec<(Kind, usize)>,  // (op_kind, length)
    pub sequence_bases: Vec<u8>,        // A, C, G, T, N
    pub num_passes: Option<i64>,        // np tag
    pub read_quality: Option<f32>,      // rq tag
    pub is_unmapped: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
    pub is_duplicate: bool,
}

impl FileHeader {
    fn from_sam_header(header: &sam::Header) -> Self {
        let mut reference_names = Vec::new();
        let mut reference_lengths = Vec::new();
        for (name, ref_seq) in header.reference_sequences() {
            reference_names.push(String::from_utf8_lossy(name.as_ref()).to_string());
            reference_lengths.push(usize::from(ref_seq.length()));
        }
        // Extract sample name from first read group if available
        let sample_name = header
            .read_groups()
            .iter()
            .next()
            .and_then(|(_, rg)| {
                rg.other_fields()
                    .get(&sam::header::record::value::map::read_group::tag::SAMPLE)
                    .map(|s| s.to_string())
            });
        FileHeader {
            reference_names,
            reference_lengths,
            sample_name,
        }
    }
}

fn extract_read_record(record: &dyn RecordFields) -> io::Result<ReadRecord> {
    // This is a helper — the actual extraction logic is inline in the
    // process_bam and process_cram functions since noodles BAM and CRAM
    // records have concrete (not trait-object) types. See the actual
    // implementation in each function below.
    unimplemented!("use format-specific extraction")
}

/// Process a BAM file, calling `callback` for each primary aligned read.
pub fn process_bam_file<F>(path: &Path, mut callback: F) -> io::Result<FileHeader>
where
    F: FnMut(ReadRecord),
{
    let mut reader = File::open(path).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;
    let file_header = FileHeader::from_sam_header(&header);

    for result in reader.records() {
        let record = result?;
        let flags = record.flags();
        let is_unmapped = flags.is_unmapped();
        let is_secondary = flags.is_secondary();
        let is_supplementary = flags.is_supplementary();
        let is_duplicate = flags.is_duplicate();

        // Skip secondary and supplementary alignments
        if is_secondary || is_supplementary {
            continue;
        }

        let mapping_quality = record
            .mapping_quality()
            .map(|mq| mq.get())
            .unwrap_or(255);

        let reference_id = record.reference_sequence_id().transpose()?;

        let alignment_start = record
            .alignment_start()
            .transpose()?
            .map(|pos| usize::from(pos) - 1); // convert 1-based to 0-based

        // Extract CIGAR ops
        let mut cigar_ops = Vec::new();
        let cigar = record.cigar();
        for op_result in cigar.iter() {
            let op = op_result?;
            cigar_ops.push((op.kind(), op.len()));
        }

        // Extract sequence bases
        let seq = record.sequence();
        let sequence_length = seq.len();
        let mut sequence_bases = Vec::with_capacity(sequence_length);
        for base in seq.iter() {
            sequence_bases.push(base);
        }

        // Extract quality scores (raw Phred, not ASCII-offset)
        let qual = record.quality_scores();
        let quality_scores: Vec<u8> = qual.iter().collect();

        // Extract HiFi tags
        let data = record.data();
        let np_tag = Tag::new(b'n', b'p');
        let num_passes: Option<i64> = data
            .get(&np_tag)
            .transpose()?
            .and_then(|v| v.as_int());

        let rq_tag = Tag::new(b'r', b'q');
        let read_quality: Option<f32> = data.get(&rq_tag).transpose()?.and_then(|v| {
            match v {
                noodles::sam::alignment::record::data::field::Value::Float(f) => Some(f),
                _ => None,
            }
        });

        callback(ReadRecord {
            sequence_length,
            quality_scores,
            mapping_quality,
            reference_id,
            alignment_start,
            cigar_ops,
            sequence_bases,
            num_passes,
            read_quality,
            is_unmapped,
            is_secondary,
            is_supplementary,
            is_duplicate,
        });
    }

    Ok(file_header)
}

/// Process a CRAM file, calling `callback` for each primary aligned read.
pub fn process_cram_file<F>(
    path: &Path,
    reference: Option<&Path>,
    mut callback: F,
) -> io::Result<FileHeader>
where
    F: FnMut(ReadRecord),
{
    let reference_repo = reference
        .map(|r| fasta::io::indexed_reader::Builder::default().build_from_path(r))
        .transpose()
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
        .map(fasta::repository::adapters::IndexedReader::new)
        .map(fasta::Repository::new)
        .unwrap_or_default();

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(reference_repo)
        .build_from_path(path)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

    let header = reader.read_header()?;
    let file_header = FileHeader::from_sam_header(&header);

    for result in reader.records(&header) {
        let record = result?;
        let flags = record.flags()?;
        let is_unmapped = flags.is_unmapped();
        let is_secondary = flags.is_secondary();
        let is_supplementary = flags.is_supplementary();
        let is_duplicate = flags.is_duplicate();

        if is_secondary || is_supplementary {
            continue;
        }

        let mapping_quality = record
            .mapping_quality()?
            .map(|mq| mq.get())
            .unwrap_or(255);

        let reference_id = record.reference_sequence_id().transpose()?;
        let alignment_start = record
            .alignment_start()
            .transpose()?
            .map(|pos| usize::from(pos) - 1);

        let mut cigar_ops = Vec::new();
        let cigar = record.cigar();
        for op_result in cigar.iter() {
            let op = op_result?;
            cigar_ops.push((op.kind(), op.len()));
        }

        let seq = record.sequence();
        let sequence_length = seq.len();
        let mut sequence_bases = Vec::with_capacity(sequence_length);
        for base in seq.iter() {
            sequence_bases.push(base);
        }

        let qual = record.quality_scores();
        let quality_scores: Vec<u8> = qual.iter().collect();

        let data = record.data();
        let np_tag = Tag::new(b'n', b'p');
        let num_passes: Option<i64> = data.get(&np_tag).transpose()?.and_then(|v| v.as_int());

        let rq_tag = Tag::new(b'r', b'q');
        let read_quality: Option<f32> = data.get(&rq_tag).transpose()?.and_then(|v| {
            match v {
                noodles::sam::alignment::record::data::field::Value::Float(f) => Some(f),
                _ => None,
            }
        });

        callback(ReadRecord {
            sequence_length,
            quality_scores,
            mapping_quality,
            reference_id,
            alignment_start,
            cigar_ops,
            sequence_bases,
            num_passes,
            read_quality,
            is_unmapped,
            is_secondary,
            is_supplementary,
            is_duplicate,
        });
    }

    Ok(file_header)
}

/// Dispatch to BAM or CRAM reader based on file extension.
pub fn process_alignment_file<F>(
    path: &Path,
    reference: Option<&Path>,
    callback: F,
) -> io::Result<FileHeader>
where
    F: FnMut(ReadRecord),
{
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();
    match ext.as_str() {
        "bam" => process_bam_file(path, callback),
        "cram" => process_cram_file(path, reference, callback),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Unsupported file format: {}", ext),
        )),
    }
}
```

**Step 2: Register the module in src/lib.rs**

Add `mod bam_reader;` to the top of `src/lib.rs`.

**Step 3: Verify it compiles**

Run: `cd /Users/arq5x/src/quinlan/hifi-qc && cargo check`
Expected: Compiles with no errors (warnings OK).

**Step 4: Commit**

```bash
git add src/bam_reader.rs src/lib.rs
git commit -m "feat: add BAM/CRAM reader module with noodles"
```

---

## Task 3: Metrics accumulator

**Files:**
- Create: `src/metrics.rs`
- Modify: `src/lib.rs`

This module defines the `SampleMetrics` struct that accumulates all QC metrics during the single pass, and a `process_read()` method that updates all accumulators from a `ReadRecord`.

**Step 1: Create src/metrics.rs**

```rust
use crate::bam_reader::ReadRecord;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::collections::HashMap;

/// Accumulated QC metrics for one sample.
pub struct SampleMetrics {
    pub sample_name: String,

    // Read counts
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub unmapped_reads: u64,
    pub duplicate_reads: u64,
    pub total_bases: u64,

    // Read length distribution (length -> count)
    pub length_histogram: HashMap<u32, u64>,
    pub lengths_for_n50: Vec<u32>, // all read lengths, for N50 computation

    // Q-score distribution (phred score -> count)
    pub qscore_histogram: HashMap<u8, u64>,

    // Read length vs quality (subsampled pairs for scatter plot)
    pub length_quality_pairs: Vec<(u32, f32)>,
    pub length_quality_subsample_count: usize,

    // HiFi passes (np tag value -> count)
    pub passes_histogram: HashMap<u32, u64>,
    // Read quality (rq tag, binned to 2 decimal places -> count)
    pub rq_histogram: HashMap<u32, u64>, // key = (rq * 10000) as u32

    // Coverage: depth arrays, one per reference sequence
    pub depth_arrays: Vec<Vec<u32>>,

    // Mapping quality distribution (MAPQ -> count)
    pub mapq_histogram: HashMap<u8, u64>,

    // GC content of reads (GC fraction binned to 1% -> count)
    pub gc_content_histogram: HashMap<u8, u64>, // key = (gc_fraction * 100) as u8

    // Mismatch analysis: substitution type -> count
    // Key format: "A>C", "A>G", "A>T", "C>A", "C>G", "C>T", etc.
    // From CIGAR X ops + sequence data (when =/X encoding used)
    pub substitution_counts: HashMap<String, u64>,
    pub total_mismatches: u64,

    // Indel analysis
    pub insertion_size_histogram: HashMap<u32, u64>, // insertion length -> count
    pub deletion_size_histogram: HashMap<u32, u64>,  // deletion length -> count
    pub total_insertions: u64,
    pub total_insertion_bases: u64,
    pub total_deletions: u64,
    pub total_deletion_bases: u64,

    // For sampling mode
    pub is_sampled: bool,
}

impl SampleMetrics {
    pub fn new(
        sample_name: String,
        reference_lengths: &[usize],
        skip_coverage: bool,
    ) -> Self {
        let depth_arrays = if skip_coverage {
            vec![]
        } else {
            reference_lengths
                .iter()
                .map(|&len| vec![0u32; len])
                .collect()
        };

        SampleMetrics {
            sample_name,
            total_reads: 0,
            mapped_reads: 0,
            unmapped_reads: 0,
            duplicate_reads: 0,
            total_bases: 0,
            length_histogram: HashMap::new(),
            lengths_for_n50: Vec::new(),
            qscore_histogram: HashMap::new(),
            length_quality_pairs: Vec::new(),
            length_quality_subsample_count: 10_000,
            passes_histogram: HashMap::new(),
            rq_histogram: HashMap::new(),
            depth_arrays,
            mapq_histogram: HashMap::new(),
            gc_content_histogram: HashMap::new(),
            substitution_counts: HashMap::new(),
            total_mismatches: 0,
            insertion_size_histogram: HashMap::new(),
            deletion_size_histogram: HashMap::new(),
            total_insertions: 0,
            total_insertion_bases: 0,
            total_deletions: 0,
            total_deletion_bases: 0,
            is_sampled: skip_coverage,
        }
    }

    /// Process a single read record, updating all accumulators.
    pub fn process_read(&mut self, read: &ReadRecord) {
        self.total_reads += 1;
        self.total_bases += read.sequence_length as u64;

        if read.is_unmapped {
            self.unmapped_reads += 1;
        } else {
            self.mapped_reads += 1;
        }
        if read.is_duplicate {
            self.duplicate_reads += 1;
        }

        // Read length
        let len = read.sequence_length as u32;
        *self.length_histogram.entry(len).or_insert(0) += 1;
        self.lengths_for_n50.push(len);

        // Mean Q-score for this read
        if !read.quality_scores.is_empty() {
            let mean_q = read.quality_scores.iter().map(|&q| q as u32).sum::<u32>()
                / read.quality_scores.len() as u32;
            *self.qscore_histogram.entry(mean_q as u8).or_insert(0) += 1;

            // Subsample length-quality pairs for scatter plot
            if self.length_quality_pairs.len() < self.length_quality_subsample_count {
                self.length_quality_pairs.push((len, mean_q as f32));
            }
        }

        // Mapping quality
        if !read.is_unmapped {
            *self.mapq_histogram.entry(read.mapping_quality).or_insert(0) += 1;
        }

        // HiFi passes (np tag)
        if let Some(np) = read.num_passes {
            *self.passes_histogram.entry(np as u32).or_insert(0) += 1;
        }

        // Read quality (rq tag)
        if let Some(rq) = read.read_quality {
            let rq_bin = (rq * 10000.0) as u32;
            *self.rq_histogram.entry(rq_bin).or_insert(0) += 1;
        }

        // GC content of the read
        if !read.sequence_bases.is_empty() {
            let gc_count = read
                .sequence_bases
                .iter()
                .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
                .count();
            let gc_pct = ((gc_count as f64 / read.sequence_bases.len() as f64) * 100.0) as u8;
            *self.gc_content_histogram.entry(gc_pct).or_insert(0) += 1;
        }

        // Skip coverage/CIGAR analysis for unmapped reads
        if read.is_unmapped {
            return;
        }

        // CIGAR walking: update depth array + count errors
        if let (Some(ref_id), Some(aln_start)) = (read.reference_id, read.alignment_start) {
            let mut ref_pos = aln_start;
            let mut _read_pos: usize = 0;

            let has_depth = !self.depth_arrays.is_empty()
                && ref_id < self.depth_arrays.len();

            for &(kind, op_len) in &read.cigar_ops {
                match kind {
                    Kind::Match | Kind::SequenceMatch => {
                        // M or =: increment depth, advance both
                        if has_depth {
                            let depth = &mut self.depth_arrays[ref_id];
                            let end = (ref_pos + op_len).min(depth.len());
                            for pos in ref_pos..end {
                                depth[pos] = depth[pos].saturating_add(1);
                            }
                        }
                        ref_pos += op_len;
                        _read_pos += op_len;
                    }
                    Kind::SequenceMismatch => {
                        // X: increment depth, count as mismatch
                        if has_depth {
                            let depth = &mut self.depth_arrays[ref_id];
                            let end = (ref_pos + op_len).min(depth.len());
                            for pos in ref_pos..end {
                                depth[pos] = depth[pos].saturating_add(1);
                            }
                        }
                        self.total_mismatches += op_len as u64;
                        ref_pos += op_len;
                        _read_pos += op_len;
                    }
                    Kind::Insertion => {
                        self.total_insertions += 1;
                        self.total_insertion_bases += op_len as u64;
                        *self
                            .insertion_size_histogram
                            .entry(op_len as u32)
                            .or_insert(0) += 1;
                        _read_pos += op_len;
                    }
                    Kind::Deletion => {
                        self.total_deletions += 1;
                        self.total_deletion_bases += op_len as u64;
                        *self
                            .deletion_size_histogram
                            .entry(op_len as u32)
                            .or_insert(0) += 1;
                        ref_pos += op_len;
                    }
                    Kind::SoftClip => {
                        _read_pos += op_len;
                    }
                    Kind::HardClip | Kind::Pad => {}
                    Kind::Skip => {
                        ref_pos += op_len;
                    }
                }
            }
        }
    }

    /// Compute N50 from collected read lengths.
    pub fn compute_n50(&self) -> u32 {
        if self.lengths_for_n50.is_empty() {
            return 0;
        }
        let mut sorted = self.lengths_for_n50.clone();
        sorted.sort_unstable_by(|a, b| b.cmp(a)); // descending
        let total: u64 = sorted.iter().map(|&l| l as u64).sum();
        let half = total / 2;
        let mut cumulative: u64 = 0;
        for &len in &sorted {
            cumulative += len as u64;
            if cumulative >= half {
                return len;
            }
        }
        sorted[0]
    }

    /// Compute mean read length.
    pub fn mean_read_length(&self) -> f64 {
        if self.total_reads == 0 {
            return 0.0;
        }
        self.total_bases as f64 / self.total_reads as f64
    }

    /// Compute median read length.
    pub fn median_read_length(&self) -> u32 {
        if self.lengths_for_n50.is_empty() {
            return 0;
        }
        let mut sorted = self.lengths_for_n50.clone();
        sorted.sort_unstable();
        let mid = sorted.len() / 2;
        if sorted.len() % 2 == 0 {
            (sorted[mid - 1] + sorted[mid]) / 2
        } else {
            sorted[mid]
        }
    }

    /// Compute mean Q-score across all reads.
    pub fn mean_qscore(&self) -> f64 {
        let total_q: u64 = self
            .qscore_histogram
            .iter()
            .map(|(&q, &count)| q as u64 * count)
            .sum();
        let total_reads: u64 = self.qscore_histogram.values().sum();
        if total_reads == 0 {
            return 0.0;
        }
        total_q as f64 / total_reads as f64
    }
}
```

**Step 2: Register in src/lib.rs**

Add `mod metrics;` to `src/lib.rs`.

**Step 3: Verify it compiles**

Run: `cd /Users/arq5x/src/quinlan/hifi-qc && cargo check`
Expected: Compiles.

**Step 4: Commit**

```bash
git add src/metrics.rs src/lib.rs
git commit -m "feat: add metrics accumulator with single-pass read processing"
```

---

## Task 4: Coverage summarization with rayon

**Files:**
- Create: `src/coverage.rs`
- Modify: `src/lib.rs`

This module takes the populated depth arrays from `SampleMetrics` and computes coverage statistics in parallel across chromosomes.

**Step 1: Create src/coverage.rs**

```rust
use rayon::prelude::*;

/// Per-chromosome coverage statistics.
pub struct ChromCoverage {
    pub name: String,
    pub length: usize,
    pub mean: f64,
    pub median: f64,
    pub stdev: f64,
}

/// Genome-wide coverage summary.
pub struct CoverageSummary {
    pub per_chrom: Vec<ChromCoverage>,
    pub genome_mean: f64,
    pub genome_median: f64,
    pub genome_stdev: f64,
    /// Coverage histogram: index = depth, value = number of bases at that depth.
    /// Capped at max_depth.
    pub coverage_histogram: Vec<u64>,
    /// Genome fraction at >= X coverage, for X in 0..max_depth.
    pub cumulative_coverage: Vec<f64>,
}

/// Compute coverage statistics from depth arrays in parallel.
pub fn summarize_coverage(
    depth_arrays: &[Vec<u32>],
    reference_names: &[String],
    max_depth: u32,
) -> CoverageSummary {
    let max_depth_usize = max_depth as usize;

    // Per-chromosome stats (parallel)
    let per_chrom: Vec<ChromCoverage> = depth_arrays
        .par_iter()
        .enumerate()
        .map(|(i, depth)| {
            let name = reference_names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("ref_{}", i));
            let n = depth.len() as f64;
            if n == 0.0 {
                return ChromCoverage {
                    name,
                    length: 0,
                    mean: 0.0,
                    median: 0.0,
                    stdev: 0.0,
                };
            }

            let sum: u64 = depth.iter().map(|&d| d as u64).sum();
            let mean = sum as f64 / n;

            // Median via histogram (avoids sorting the full array)
            let mut hist = vec![0u64; max_depth_usize + 1];
            for &d in depth {
                let bin = (d as usize).min(max_depth_usize);
                hist[bin] += 1;
            }
            let median = histogram_median(&hist, depth.len());

            // Stdev
            let variance: f64 = depth
                .iter()
                .map(|&d| {
                    let diff = d as f64 - mean;
                    diff * diff
                })
                .sum::<f64>()
                / n;
            let stdev = variance.sqrt();

            ChromCoverage {
                name,
                length: depth.len(),
                mean,
                median,
                stdev,
            }
        })
        .collect();

    // Genome-wide coverage histogram (parallel merge)
    let coverage_histogram: Vec<u64> = depth_arrays
        .par_iter()
        .map(|depth| {
            let mut hist = vec![0u64; max_depth_usize + 1];
            for &d in depth {
                let bin = (d as usize).min(max_depth_usize);
                hist[bin] += 1;
            }
            hist
        })
        .reduce(
            || vec![0u64; max_depth_usize + 1],
            |mut a, b| {
                for (i, v) in b.iter().enumerate() {
                    a[i] += v;
                }
                a
            },
        );

    let total_bases: u64 = coverage_histogram.iter().sum();

    // Genome-wide mean
    let genome_mean = if total_bases > 0 {
        coverage_histogram
            .iter()
            .enumerate()
            .map(|(depth, &count)| depth as f64 * count as f64)
            .sum::<f64>()
            / total_bases as f64
    } else {
        0.0
    };

    // Genome-wide median from histogram
    let genome_median = histogram_median(&coverage_histogram, total_bases as usize);

    // Genome-wide stdev from histogram
    let genome_stdev = if total_bases > 0 {
        let variance: f64 = coverage_histogram
            .iter()
            .enumerate()
            .map(|(depth, &count)| {
                let diff = depth as f64 - genome_mean;
                diff * diff * count as f64
            })
            .sum::<f64>()
            / total_bases as f64;
        variance.sqrt()
    } else {
        0.0
    };

    // Cumulative coverage: fraction of genome at >= X depth
    let cumulative_coverage: Vec<f64> = (0..=max_depth_usize)
        .map(|threshold| {
            if total_bases == 0 {
                return 0.0;
            }
            let bases_at_or_above: u64 = coverage_histogram[threshold..].iter().sum();
            bases_at_or_above as f64 / total_bases as f64
        })
        .collect();

    CoverageSummary {
        per_chrom,
        genome_mean,
        genome_median,
        genome_stdev,
        coverage_histogram,
        cumulative_coverage,
    }
}

/// Compute median from a histogram of counts.
fn histogram_median(hist: &[u64], total: usize) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let half = total / 2;
    let mut cumulative: u64 = 0;
    for (depth, &count) in hist.iter().enumerate() {
        cumulative += count;
        if cumulative as usize > half {
            return depth as f64;
        }
    }
    0.0
}
```

**Step 2: Register in src/lib.rs**

Add `mod coverage;` to `src/lib.rs`.

**Step 3: Verify it compiles**

Run: `cd /Users/arq5x/src/quinlan/hifi-qc && cargo check`
Expected: Compiles.

**Step 4: Commit**

```bash
git add src/coverage.rs src/lib.rs
git commit -m "feat: add parallel coverage summarization with rayon"
```

---

## Task 5: GC bias computation

**Files:**
- Create: `src/gc_bias.rs`
- Modify: `src/lib.rs`

Computes reference GC content in windows, then stratifies observed coverage by GC bin.

**Step 1: Create src/gc_bias.rs**

```rust
use noodles::fasta;
use rayon::prelude::*;
use std::io;
use std::path::Path;

const WINDOW_SIZE: usize = 100;
const N_GC_BINS: usize = 101; // 0% to 100% inclusive

/// GC bias result: expected vs. observed coverage by GC bin.
pub struct GcBiasResult {
    /// For each GC bin (0..=100): number of reference windows in that bin.
    pub reference_gc_distribution: Vec<u64>,
    /// For each GC bin (0..=100): mean coverage of windows in that bin.
    pub coverage_by_gc: Vec<f64>,
    /// The observed GC content distribution of reads (from SampleMetrics).
    pub read_gc_distribution: Vec<u64>,
}

/// Compute GC bias by stratifying coverage by reference GC content.
///
/// - `reference_path`: path to the reference FASTA (must be indexed with .fai)
/// - `depth_arrays`: per-chromosome depth arrays from coverage computation
/// - `reference_names`: chromosome names in BAM header order
pub fn compute_gc_bias(
    reference_path: &Path,
    depth_arrays: &[Vec<u32>],
    reference_names: &[String],
) -> io::Result<GcBiasResult> {
    // Read reference sequences
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference_path)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

    // For each chromosome, compute per-window GC content and mean coverage
    // Collect (gc_bin, mean_coverage) pairs for all windows across all chromosomes
    let mut gc_coverage_sums = vec![0.0f64; N_GC_BINS];
    let mut gc_window_counts = vec![0u64; N_GC_BINS];

    for (chrom_idx, chrom_name) in reference_names.iter().enumerate() {
        if chrom_idx >= depth_arrays.len() {
            continue;
        }
        let depth = &depth_arrays[chrom_idx];
        if depth.is_empty() {
            continue;
        }

        // Read the reference sequence for this chromosome
        let region_str = format!("{}", chrom_name);
        let region: noodles::core::Region = region_str
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("{}", e)))?;

        let record = reader
            .query(&region)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

        let seq: Vec<u8> = record
            .sequence()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();

        // Process in windows
        let n_windows = depth.len() / WINDOW_SIZE;
        for w in 0..n_windows {
            let start = w * WINDOW_SIZE;
            let end = start + WINDOW_SIZE;

            // GC content of this window
            let gc_count = seq[start..end]
                .iter()
                .filter(|&&b| b == b'G' || b == b'C')
                .count();
            let n_count = seq[start..end]
                .iter()
                .filter(|&&b| b == b'N')
                .count();

            // Skip windows with too many Ns
            if n_count > WINDOW_SIZE / 2 {
                continue;
            }

            let gc_pct = ((gc_count as f64 / (WINDOW_SIZE - n_count) as f64) * 100.0) as usize;
            let gc_bin = gc_pct.min(100);

            // Mean coverage in this window
            let window_cov: f64 = depth[start..end].iter().map(|&d| d as f64).sum::<f64>()
                / WINDOW_SIZE as f64;

            gc_coverage_sums[gc_bin] += window_cov;
            gc_window_counts[gc_bin] += 1;
        }
    }

    // Compute mean coverage per GC bin
    let coverage_by_gc: Vec<f64> = gc_coverage_sums
        .iter()
        .zip(gc_window_counts.iter())
        .map(|(&sum, &count)| if count > 0 { sum / count as f64 } else { 0.0 })
        .collect();

    Ok(GcBiasResult {
        reference_gc_distribution: gc_window_counts,
        coverage_by_gc,
        read_gc_distribution: vec![], // filled in by caller from SampleMetrics
    })
}
```

**Step 2: Register in src/lib.rs**

Add `mod gc_bias;` to `src/lib.rs`.

**Step 3: Verify it compiles**

Run: `cd /Users/arq5x/src/quinlan/hifi-qc && cargo check`
Expected: Compiles.

**Step 4: Commit**

```bash
git add src/gc_bias.rs src/lib.rs
git commit -m "feat: add GC bias computation from reference FASTA"
```

---

## Task 6: Wire up PyO3 entry point

**Files:**
- Modify: `src/lib.rs`

Connect all Rust modules into the `process_bam_files` PyO3 function. This function:
1. Sets up rayon thread pool
2. For each BAM/CRAM: single pass collecting all metrics
3. Summarizes coverage (parallel)
4. Computes GC bias (if reference provided)
5. Returns results as Python dicts

**Step 1: Rewrite src/lib.rs**

```rust
mod bam_reader;
mod coverage;
mod gc_bias;
mod metrics;

use bam_reader::{process_alignment_file, ReadRecord};
use coverage::summarize_coverage;
use gc_bias::compute_gc_bias;
use metrics::SampleMetrics;

use pyo3::prelude::*;
use pyo3::types::PyDict;
use rand::Rng;
use rand::SeedableRng;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::path::Path;

const MAX_COVERAGE_DEPTH: u32 = 1000;

#[pyfunction]
#[pyo3(signature = (bam_paths, reference=None, threads=None, sample_fraction=None, seed=None, sample_names=None))]
fn process_bam_files(
    py: Python<'_>,
    bam_paths: Vec<String>,
    reference: Option<String>,
    threads: Option<usize>,
    sample_fraction: Option<f64>,
    seed: Option<u64>,
    sample_names: Option<Vec<String>>,
) -> PyResult<PyObject> {
    // Set up thread pool
    if let Some(n) = threads {
        ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok(); // ignore error if already built
    }

    let ref_path = reference.as_deref().map(Path::new);
    let is_sampling = sample_fraction.map_or(false, |f| f < 1.0);

    let mut all_results: Vec<PyObject> = Vec::new();

    for (file_idx, bam_path) in bam_paths.iter().enumerate() {
        let path = Path::new(bam_path);

        // Set up optional RNG for sampling
        let mut rng = if is_sampling {
            let s = seed.unwrap_or(42) + file_idx as u64;
            Some(rand::rngs::StdRng::seed_from_u64(s))
        } else {
            None
        };
        let fraction = sample_fraction.unwrap_or(1.0);

        // Create metrics accumulator (skip coverage if sampling)
        let header = {
            // We need the header first to allocate depth arrays
            // Process the file and collect metrics
            let mut sample_metrics: Option<SampleMetrics> = None;

            let file_header = process_alignment_file(path, ref_path, |read: ReadRecord| {
                // Sampling: probabilistically skip reads
                if let Some(ref mut r) = rng {
                    if r.random::<f64>() >= fraction {
                        return;
                    }
                }

                // Lazy-init metrics on first read (we need header info first)
                // This is handled by the closure capturing sample_metrics
                if let Some(ref mut m) = sample_metrics {
                    m.process_read(&read);
                }
            })
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

            // If sample_metrics was never initialized, create it now
            // (This happens because we need to init it before the first callback)
            // We need to restructure: init metrics before iterating.
            file_header
        };

        // Actually, we need to restructure: read header first, then iterate.
        // Let's do a two-phase approach within the same function.
        // Phase 1: open file, read header
        // Phase 2: iterate records with metrics accumulator

        // Re-do with proper structure:
        let sample_name = sample_names
            .as_ref()
            .and_then(|names| names.get(file_idx).cloned())
            .or(header.sample_name.clone())
            .unwrap_or_else(|| {
                Path::new(bam_path)
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .to_string()
            });

        let mut m = SampleMetrics::new(
            sample_name,
            &header.reference_lengths,
            is_sampling,
        );

        // Set up RNG again for the actual pass
        let mut rng = if is_sampling {
            let s = seed.unwrap_or(42) + file_idx as u64;
            Some(rand::rngs::StdRng::seed_from_u64(s))
        } else {
            None
        };

        // Process file (this is the actual single pass)
        process_alignment_file(path, ref_path, |read: ReadRecord| {
            if let Some(ref mut r) = rng {
                if r.random::<f64>() >= fraction {
                    return;
                }
            }
            m.process_read(&read);
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        // Summarize coverage (parallel, only if not sampling)
        let coverage = if !is_sampling && !m.depth_arrays.is_empty() {
            Some(summarize_coverage(
                &m.depth_arrays,
                &header.reference_names,
                MAX_COVERAGE_DEPTH,
            ))
        } else {
            None
        };

        // Compute GC bias (if reference provided and not sampling)
        let gc_bias = if let (Some(ref_p), false) = (ref_path, is_sampling) {
            compute_gc_bias(ref_p, &m.depth_arrays, &header.reference_names).ok()
        } else {
            None
        };

        // Build Python dict for this sample
        let result = build_sample_dict(py, &m, &header, coverage.as_ref(), gc_bias.as_ref())?;
        all_results.push(result);
    }

    Ok(all_results.into_pyobject(py)?.into_any().unbind())
}

fn build_sample_dict(
    py: Python<'_>,
    m: &SampleMetrics,
    header: &bam_reader::FileHeader,
    coverage: Option<&coverage::CoverageSummary>,
    gc_bias: Option<&gc_bias::GcBiasResult>,
) -> PyResult<PyObject> {
    let dict = PyDict::new(py);

    // Summary stats
    dict.set_item("sample_name", &m.sample_name)?;
    dict.set_item("total_reads", m.total_reads)?;
    dict.set_item("mapped_reads", m.mapped_reads)?;
    dict.set_item("unmapped_reads", m.unmapped_reads)?;
    dict.set_item("duplicate_reads", m.duplicate_reads)?;
    dict.set_item("total_bases", m.total_bases)?;
    dict.set_item("mean_read_length", m.mean_read_length())?;
    dict.set_item("median_read_length", m.median_read_length())?;
    dict.set_item("n50_read_length", m.compute_n50())?;
    dict.set_item("mean_qscore", m.mean_qscore())?;
    dict.set_item("is_sampled", m.is_sampled)?;

    // Read length histogram (sorted by length)
    let mut length_hist: Vec<(u32, u64)> = m.length_histogram.iter().map(|(&k, &v)| (k, v)).collect();
    length_hist.sort_by_key(|&(k, _)| k);
    let length_keys: Vec<u32> = length_hist.iter().map(|&(k, _)| k).collect();
    let length_vals: Vec<u64> = length_hist.iter().map(|&(_, v)| v).collect();
    dict.set_item("length_hist_keys", length_keys)?;
    dict.set_item("length_hist_values", length_vals)?;

    // Q-score histogram
    let mut qscore_hist: Vec<(u8, u64)> = m.qscore_histogram.iter().map(|(&k, &v)| (k, v)).collect();
    qscore_hist.sort_by_key(|&(k, _)| k);
    let q_keys: Vec<u8> = qscore_hist.iter().map(|&(k, _)| k).collect();
    let q_vals: Vec<u64> = qscore_hist.iter().map(|&(_, v)| v).collect();
    dict.set_item("qscore_hist_keys", q_keys)?;
    dict.set_item("qscore_hist_values", q_vals)?;

    // Length vs quality pairs
    let lq_lengths: Vec<u32> = m.length_quality_pairs.iter().map(|&(l, _)| l).collect();
    let lq_quals: Vec<f32> = m.length_quality_pairs.iter().map(|&(_, q)| q).collect();
    dict.set_item("length_quality_lengths", lq_lengths)?;
    dict.set_item("length_quality_scores", lq_quals)?;

    // HiFi passes histogram
    let mut passes_hist: Vec<(u32, u64)> = m.passes_histogram.iter().map(|(&k, &v)| (k, v)).collect();
    passes_hist.sort_by_key(|&(k, _)| k);
    let p_keys: Vec<u32> = passes_hist.iter().map(|&(k, _)| k).collect();
    let p_vals: Vec<u64> = passes_hist.iter().map(|&(_, v)| v).collect();
    dict.set_item("passes_hist_keys", p_keys)?;
    dict.set_item("passes_hist_values", p_vals)?;

    // Read quality (rq) histogram
    let mut rq_hist: Vec<(u32, u64)> = m.rq_histogram.iter().map(|(&k, &v)| (k, v)).collect();
    rq_hist.sort_by_key(|&(k, _)| k);
    let rq_keys: Vec<f64> = rq_hist.iter().map(|&(k, _)| k as f64 / 10000.0).collect();
    let rq_vals: Vec<u64> = rq_hist.iter().map(|&(_, v)| v).collect();
    dict.set_item("rq_hist_keys", rq_keys)?;
    dict.set_item("rq_hist_values", rq_vals)?;

    // Mapping quality histogram
    let mut mapq_hist: Vec<(u8, u64)> = m.mapq_histogram.iter().map(|(&k, &v)| (k, v)).collect();
    mapq_hist.sort_by_key(|&(k, _)| k);
    let mq_keys: Vec<u8> = mapq_hist.iter().map(|&(k, _)| k).collect();
    let mq_vals: Vec<u64> = mapq_hist.iter().map(|&(_, v)| v).collect();
    dict.set_item("mapq_hist_keys", mq_keys)?;
    dict.set_item("mapq_hist_values", mq_vals)?;

    // GC content histogram
    let mut gc_hist: Vec<(u8, u64)> = m.gc_content_histogram.iter().map(|(&k, &v)| (k, v)).collect();
    gc_hist.sort_by_key(|&(k, _)| k);
    let gc_keys: Vec<u8> = gc_hist.iter().map(|&(k, _)| k).collect();
    let gc_vals: Vec<u64> = gc_hist.iter().map(|&(_, v)| v).collect();
    dict.set_item("gc_content_keys", gc_keys)?;
    dict.set_item("gc_content_values", gc_vals)?;

    // Mismatch / indel summary
    dict.set_item("total_mismatches", m.total_mismatches)?;
    dict.set_item("total_insertions", m.total_insertions)?;
    dict.set_item("total_insertion_bases", m.total_insertion_bases)?;
    dict.set_item("total_deletions", m.total_deletions)?;
    dict.set_item("total_deletion_bases", m.total_deletion_bases)?;

    // Substitution counts
    let sub_keys: Vec<String> = m.substitution_counts.keys().cloned().collect();
    let sub_vals: Vec<u64> = sub_keys.iter().map(|k| m.substitution_counts[k]).collect();
    dict.set_item("substitution_types", sub_keys)?;
    dict.set_item("substitution_counts", sub_vals)?;

    // Insertion size histogram
    let mut ins_hist: Vec<(u32, u64)> = m.insertion_size_histogram.iter().map(|(&k, &v)| (k, v)).collect();
    ins_hist.sort_by_key(|&(k, _)| k);
    let ins_keys: Vec<u32> = ins_hist.iter().map(|&(k, _)| k).collect();
    let ins_vals: Vec<u64> = ins_hist.iter().map(|&(_, v)| v).collect();
    dict.set_item("insertion_size_keys", ins_keys)?;
    dict.set_item("insertion_size_values", ins_vals)?;

    // Deletion size histogram
    let mut del_hist: Vec<(u32, u64)> = m.deletion_size_histogram.iter().map(|(&k, &v)| (k, v)).collect();
    del_hist.sort_by_key(|&(k, _)| k);
    let del_keys: Vec<u32> = del_hist.iter().map(|&(k, _)| k).collect();
    let del_vals: Vec<u64> = del_hist.iter().map(|&(_, v)| v).collect();
    dict.set_item("deletion_size_keys", del_keys)?;
    dict.set_item("deletion_size_values", del_vals)?;

    // Coverage data (only if computed)
    if let Some(cov) = coverage {
        dict.set_item("genome_mean_coverage", cov.genome_mean)?;
        dict.set_item("genome_median_coverage", cov.genome_median)?;
        dict.set_item("genome_stdev_coverage", cov.genome_stdev)?;

        let chrom_names: Vec<&str> = cov.per_chrom.iter().map(|c| c.name.as_str()).collect();
        let chrom_means: Vec<f64> = cov.per_chrom.iter().map(|c| c.mean).collect();
        let chrom_lengths: Vec<usize> = cov.per_chrom.iter().map(|c| c.length).collect();
        dict.set_item("chrom_names", chrom_names)?;
        dict.set_item("chrom_mean_coverages", chrom_means)?;
        dict.set_item("chrom_lengths", chrom_lengths)?;

        dict.set_item("coverage_histogram", &cov.coverage_histogram)?;
        dict.set_item("cumulative_coverage", &cov.cumulative_coverage)?;
    }

    // GC bias data (only if computed)
    if let Some(gc) = gc_bias {
        dict.set_item("gc_bias_ref_distribution", &gc.reference_gc_distribution)?;
        dict.set_item("gc_bias_coverage_by_gc", &gc.coverage_by_gc)?;
    }

    Ok(dict.into())
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(process_bam_files, m)?)?;
    Ok(())
}
```

**Note:** The above has a structural issue — we iterate the file twice (once for header, once for data). This needs to be restructured so that `process_alignment_file` first returns the header, then we create the metrics accumulator, then iterate. Modify `bam_reader.rs` to split into a `read_header` function and a `process_records` function, OR restructure `process_alignment_file` to accept a factory closure that creates the metrics after reading the header. The simplest fix: change `process_alignment_file` to return both the header and call the callback, by reading the header first internally. The current design already does this — the callback is called during iteration which happens after the header is read. The issue is that `SampleMetrics` needs header info to allocate depth arrays, but the callback captures `&mut SampleMetrics`. Solution: split into `open_alignment_file()` -> `(header, RecordIterator)` and iterate separately. Or: pass reference_lengths into `process_alignment_file` as a separate step.

**Better approach:** Refactor `bam_reader.rs` to have two functions:
1. `read_header(path, reference)` -> `FileHeader`
2. `process_records(path, reference, callback)` -> calls callback for each record

Then in `lib.rs`:
```rust
let header = bam_reader::read_header(path, ref_path)?;
let mut m = SampleMetrics::new(sample_name, &header.reference_lengths, is_sampling);
bam_reader::process_records(path, ref_path, |read| {
    // sampling logic
    m.process_read(&read);
})?;
```

This does open the file twice but that's negligible overhead. The actual BAM iteration (the expensive part) only happens once.

**Step 2: Update bam_reader.rs to add read_header function**

Add a `read_header` function that opens the file, reads just the header, and returns `FileHeader`:

```rust
pub fn read_header(path: &Path, reference: Option<&Path>) -> io::Result<FileHeader> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_lowercase();
    match ext.as_str() {
        "bam" => {
            let mut reader = File::open(path).map(bam::io::Reader::new)?;
            let header = reader.read_header()?;
            Ok(FileHeader::from_sam_header(&header))
        }
        "cram" => {
            let reference_repo = reference
                .map(|r| fasta::io::indexed_reader::Builder::default().build_from_path(r))
                .transpose()
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
                .map(fasta::repository::adapters::IndexedReader::new)
                .map(fasta::Repository::new)
                .unwrap_or_default();
            let mut reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(reference_repo)
                .build_from_path(path)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            let header = reader.read_header()?;
            Ok(FileHeader::from_sam_header(&header))
        }
        _ => Err(io::Error::new(io::ErrorKind::InvalidInput, format!("Unsupported: {}", ext))),
    }
}
```

**Step 3: Verify it compiles**

Run: `cd /Users/arq5x/src/quinlan/hifi-qc && cargo check`
Expected: Compiles.

**Step 4: Build with maturin and test end-to-end**

Run: `cd /Users/arq5x/src/quinlan/hifi-qc && maturin develop`
Expected: Builds successfully.

Run: `hifi-qc --help`
Expected: Shows all CLI options.

**Step 5: Commit**

```bash
git add src/lib.rs src/bam_reader.rs
git commit -m "feat: wire up PyO3 entry point connecting all Rust modules"
```

---

## Task 7: HTML report template

**Files:**
- Create: `python/hifi_qc/templates/report.html.j2`
- Create: `python/hifi_qc/report.py`
- Modify: `python/hifi_qc/cli.py`

This creates the Jinja2 template and Python report generator. The template is a self-contained HTML file with Plotly.js inlined, producing a MultiQC-style interactive report.

**Step 1: Create python/hifi_qc/report.py**

```python
"""HTML report generator for hifi-qc."""

import json
from pathlib import Path

from jinja2 import Environment, FileSystemLoader


# MultiQC-style color palette
COLORS = [
    "#e41a1c",  # red
    "#377eb8",  # blue
    "#4daf4a",  # green
    "#984ea3",  # purple
    "#ff7f00",  # orange
    "#a65628",  # brown
    "#f781bf",  # pink
    "#999999",  # gray
]


def generate_report(samples: list[dict], output_path: str) -> None:
    """Generate a self-contained HTML report from sample metrics.

    Args:
        samples: List of dicts, one per sample, as returned by the Rust core.
        output_path: Path to write the HTML report.
    """
    template_dir = Path(__file__).parent / "templates"
    env = Environment(loader=FileSystemLoader(str(template_dir)), autoescape=False)
    template = env.get_template("report.html.j2")

    # Assign colors to samples
    for i, sample in enumerate(samples):
        sample["color"] = COLORS[i % len(COLORS)]

    html = template.render(
        samples=samples,
        samples_json=json.dumps(samples, default=_json_serializer),
        colors=COLORS,
        has_coverage=any("genome_mean_coverage" in s for s in samples),
        has_gc_bias=any("gc_bias_coverage_by_gc" in s for s in samples),
    )

    Path(output_path).write_text(html)


def _json_serializer(obj):
    """Handle types that json.dumps can't serialize."""
    if hasattr(obj, "__float__"):
        return float(obj)
    if hasattr(obj, "__int__"):
        return int(obj)
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")
```

**Step 2: Create python/hifi_qc/templates/report.html.j2**

This is a large template. Key sections:
- HTML head with inlined Plotly.js (loaded from CDN via script tag — the report is "self-contained" in that it has no custom external dependencies, but Plotly.js is loaded from CDN for size reasons, with a fallback message)
- Navigation sidebar
- General Statistics table
- Each QC section with a Plotly chart
- JavaScript to render all charts from the embedded JSON data

```html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>hifi-qc Report</title>
    <script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            color: #333;
            background: #f8f9fa;
            display: flex;
        }
        /* Sidebar */
        .sidebar {
            width: 220px;
            background: #2c3e50;
            color: #ecf0f1;
            position: fixed;
            top: 0;
            left: 0;
            height: 100vh;
            overflow-y: auto;
            padding: 20px 0;
            z-index: 100;
        }
        .sidebar h1 {
            font-size: 1.3em;
            padding: 0 20px 15px;
            border-bottom: 1px solid #34495e;
            margin-bottom: 10px;
        }
        .sidebar a {
            display: block;
            color: #bdc3c7;
            text-decoration: none;
            padding: 8px 20px;
            font-size: 0.85em;
            transition: background 0.2s;
        }
        .sidebar a:hover { background: #34495e; color: #fff; }
        .sidebar a.active { background: #3498db; color: #fff; }
        /* Main content */
        .main {
            margin-left: 220px;
            padding: 30px;
            width: calc(100% - 220px);
        }
        .header {
            background: #fff;
            border-radius: 8px;
            padding: 25px 30px;
            margin-bottom: 25px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        .header h1 { font-size: 1.8em; margin-bottom: 5px; }
        .header p { color: #666; }
        /* Sections */
        .section {
            background: #fff;
            border-radius: 8px;
            padding: 25px 30px;
            margin-bottom: 25px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        .section h2 {
            font-size: 1.3em;
            margin-bottom: 5px;
            color: #2c3e50;
        }
        .section .description {
            color: #666;
            font-size: 0.9em;
            margin-bottom: 15px;
        }
        .plot { width: 100%; min-height: 400px; }
        /* Summary table */
        .summary-table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9em;
        }
        .summary-table th {
            background: #2c3e50;
            color: #fff;
            padding: 10px 12px;
            text-align: left;
            font-weight: 500;
        }
        .summary-table td {
            padding: 8px 12px;
            border-bottom: 1px solid #eee;
        }
        .summary-table tr:hover { background: #f5f6fa; }
        .sample-dot {
            display: inline-block;
            width: 10px;
            height: 10px;
            border-radius: 50%;
            margin-right: 6px;
        }
        .badge {
            display: inline-block;
            background: #e74c3c;
            color: #fff;
            font-size: 0.75em;
            padding: 2px 6px;
            border-radius: 3px;
            margin-left: 6px;
        }
        .badge.sampling { background: #f39c12; }
    </style>
</head>
<body>
    <nav class="sidebar">
        <h1>hifi-qc</h1>
        <a href="#summary">General Statistics</a>
        <a href="#read-length">Read Length</a>
        <a href="#qscore">Q-Score</a>
        <a href="#length-vs-quality">Length vs Quality</a>
        <a href="#hifi-passes">HiFi Passes</a>
        {% if has_coverage %}
        <a href="#coverage-genome">Genome Coverage</a>
        <a href="#coverage-chrom">Per-Chromosome</a>
        <a href="#coverage-uniformity">Coverage Uniformity</a>
        {% endif %}
        {% if has_gc_bias %}
        <a href="#gc-bias">GC Bias</a>
        {% endif %}
        <a href="#mapq">Mapping Quality</a>
        <a href="#mismatches">Mismatches &amp; Indels</a>
    </nav>

    <div class="main">
        <div class="header">
            <h1>hifi-qc Report</h1>
            <p>{{ samples | length }} sample(s) analyzed
            {% if samples[0].is_sampled %}<span class="badge sampling">Sampling mode</span>{% endif %}
            </p>
        </div>

        <!-- General Statistics -->
        <div class="section" id="summary">
            <h2>General Statistics</h2>
            <table class="summary-table">
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>Total Reads</th>
                        <th>Total Bases</th>
                        <th>Mean Length</th>
                        <th>N50</th>
                        <th>Mean Q</th>
                        {% if has_coverage %}<th>Mean Cov</th>{% endif %}
                        <th>Mapped %</th>
                    </tr>
                </thead>
                <tbody>
                {% for s in samples %}
                    <tr>
                        <td><span class="sample-dot" style="background:{{ s.color }}"></span>{{ s.sample_name }}</td>
                        <td>{{ "{:,}".format(s.total_reads) }}</td>
                        <td>{{ "{:,.0f}".format(s.total_bases / 1e9) }} Gb</td>
                        <td>{{ "{:,.0f}".format(s.mean_read_length) }}</td>
                        <td>{{ "{:,}".format(s.n50_read_length) }}</td>
                        <td>Q{{ "{:.1f}".format(s.mean_qscore) }}</td>
                        {% if has_coverage %}<td>{{ "{:.1f}".format(s.get("genome_mean_coverage", 0)) }}x</td>{% endif %}
                        <td>{{ "{:.1f}".format(s.mapped_reads / s.total_reads * 100 if s.total_reads > 0 else 0) }}%</td>
                    </tr>
                {% endfor %}
                </tbody>
            </table>
        </div>

        <!-- Read Length Distribution -->
        <div class="section" id="read-length">
            <h2>Read Length Distribution</h2>
            <p class="description">Distribution of HiFi read lengths across all samples.</p>
            <div id="plot-read-length" class="plot"></div>
        </div>

        <!-- Q-Score Distribution -->
        <div class="section" id="qscore">
            <h2>Q-Score Distribution</h2>
            <p class="description">Distribution of per-read mean quality scores (Phred scale).</p>
            <div id="plot-qscore" class="plot"></div>
        </div>

        <!-- Read Length vs Quality -->
        <div class="section" id="length-vs-quality">
            <h2>Read Length vs Quality</h2>
            <p class="description">Relationship between read length and mean quality score.</p>
            <div id="plot-length-quality" class="plot"></div>
        </div>

        <!-- HiFi Passes -->
        <div class="section" id="hifi-passes">
            <h2>HiFi Passes</h2>
            <p class="description">Number of subreads passes per CCS molecule (np tag).</p>
            <div id="plot-passes" class="plot"></div>
        </div>

        {% if has_coverage %}
        <!-- Genome-wide Coverage -->
        <div class="section" id="coverage-genome">
            <h2>Genome-wide Coverage Distribution</h2>
            <p class="description">Distribution of per-base sequencing depth across the genome.</p>
            <div id="plot-coverage-genome" class="plot"></div>
        </div>

        <!-- Per-Chromosome Coverage -->
        <div class="section" id="coverage-chrom">
            <h2>Per-Chromosome Coverage</h2>
            <p class="description">Mean coverage depth for each reference sequence.</p>
            <div id="plot-coverage-chrom" class="plot"></div>
        </div>

        <!-- Coverage Uniformity -->
        <div class="section" id="coverage-uniformity">
            <h2>Coverage Uniformity</h2>
            <p class="description">Fraction of genome covered at or above each depth threshold.</p>
            <div id="plot-coverage-uniformity" class="plot"></div>
        </div>
        {% endif %}

        {% if has_gc_bias %}
        <!-- GC Bias -->
        <div class="section" id="gc-bias">
            <h2>GC Bias</h2>
            <p class="description">Normalized coverage as a function of GC content.</p>
            <div id="plot-gc-bias" class="plot"></div>
        </div>
        {% endif %}

        <!-- Mapping Quality -->
        <div class="section" id="mapq">
            <h2>Mapping Quality Distribution</h2>
            <p class="description">Distribution of mapping quality scores (MAPQ).</p>
            <div id="plot-mapq" class="plot"></div>
        </div>

        <!-- Mismatches & Indels -->
        <div class="section" id="mismatches">
            <h2>Mismatches &amp; Indels</h2>
            <p class="description">Insertion and deletion size distributions from CIGAR analysis.</p>
            <div id="plot-indels" class="plot"></div>
        </div>
    </div>

    <script>
    const DATA = {{ samples_json }};
    const plotConfig = {responsive: true, displayModeBar: true, displaylogo: false};
    const defaultLayout = {
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: {family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif', size: 12},
        margin: {l: 60, r: 30, t: 30, b: 50},
        legend: {orientation: 'h', y: -0.15},
        xaxis: {gridcolor: '#eee'},
        yaxis: {gridcolor: '#eee'},
    };

    // Read Length Distribution
    (function() {
        const traces = DATA.map(s => ({
            x: s.length_hist_keys,
            y: s.length_hist_values,
            type: 'bar',
            name: s.sample_name,
            marker: {color: s.color},
            opacity: 0.7,
        }));
        Plotly.newPlot('plot-read-length', traces, {
            ...defaultLayout,
            barmode: 'overlay',
            xaxis: {...defaultLayout.xaxis, title: 'Read Length (bp)'},
            yaxis: {...defaultLayout.yaxis, title: 'Count'},
        }, plotConfig);
    })();

    // Q-Score Distribution
    (function() {
        const traces = DATA.map(s => ({
            x: s.qscore_hist_keys,
            y: s.qscore_hist_values,
            type: 'bar',
            name: s.sample_name,
            marker: {color: s.color},
            opacity: 0.7,
        }));
        Plotly.newPlot('plot-qscore', traces, {
            ...defaultLayout,
            barmode: 'overlay',
            xaxis: {...defaultLayout.xaxis, title: 'Mean Q-Score (Phred)'},
            yaxis: {...defaultLayout.yaxis, title: 'Count'},
        }, plotConfig);
    })();

    // Length vs Quality scatter
    (function() {
        const traces = DATA.map(s => ({
            x: s.length_quality_lengths,
            y: s.length_quality_scores,
            mode: 'markers',
            type: 'scattergl',
            name: s.sample_name,
            marker: {color: s.color, size: 3, opacity: 0.3},
        }));
        Plotly.newPlot('plot-length-quality', traces, {
            ...defaultLayout,
            xaxis: {...defaultLayout.xaxis, title: 'Read Length (bp)'},
            yaxis: {...defaultLayout.yaxis, title: 'Mean Q-Score'},
        }, plotConfig);
    })();

    // HiFi Passes
    (function() {
        const traces = DATA.map(s => ({
            x: s.passes_hist_keys,
            y: s.passes_hist_values,
            type: 'bar',
            name: s.sample_name,
            marker: {color: s.color},
            opacity: 0.7,
        }));
        Plotly.newPlot('plot-passes', traces, {
            ...defaultLayout,
            barmode: 'overlay',
            xaxis: {...defaultLayout.xaxis, title: 'Number of Passes'},
            yaxis: {...defaultLayout.yaxis, title: 'Count'},
        }, plotConfig);
    })();

    // Coverage sections (only if data exists)
    {% if has_coverage %}
    // Genome-wide coverage distribution
    (function() {
        const traces = DATA.filter(s => s.coverage_histogram).map(s => {
            // Trim trailing zeros for cleaner plot
            let hist = s.coverage_histogram;
            let maxIdx = hist.length - 1;
            while (maxIdx > 0 && hist[maxIdx] === 0) maxIdx--;
            const trimmed = hist.slice(0, Math.min(maxIdx + 10, hist.length));
            const x = trimmed.map((_, i) => i);
            return {
                x: x,
                y: trimmed,
                type: 'bar',
                name: s.sample_name,
                marker: {color: s.color},
                opacity: 0.7,
            };
        });
        Plotly.newPlot('plot-coverage-genome', traces, {
            ...defaultLayout,
            barmode: 'overlay',
            xaxis: {...defaultLayout.xaxis, title: 'Coverage Depth', range: [0, 100]},
            yaxis: {...defaultLayout.yaxis, title: 'Number of Bases'},
        }, plotConfig);
    })();

    // Per-chromosome coverage
    (function() {
        const traces = DATA.filter(s => s.chrom_names).map(s => ({
            x: s.chrom_names,
            y: s.chrom_mean_coverages,
            type: 'bar',
            name: s.sample_name,
            marker: {color: s.color},
        }));
        Plotly.newPlot('plot-coverage-chrom', traces, {
            ...defaultLayout,
            barmode: 'group',
            xaxis: {...defaultLayout.xaxis, title: 'Chromosome', tickangle: -45},
            yaxis: {...defaultLayout.yaxis, title: 'Mean Coverage'},
        }, plotConfig);
    })();

    // Coverage uniformity
    (function() {
        const traces = DATA.filter(s => s.cumulative_coverage).map(s => {
            const x = s.cumulative_coverage.map((_, i) => i);
            return {
                x: x,
                y: s.cumulative_coverage.map(v => v * 100),
                type: 'scatter',
                mode: 'lines',
                name: s.sample_name,
                line: {color: s.color, width: 2},
            };
        });
        Plotly.newPlot('plot-coverage-uniformity', traces, {
            ...defaultLayout,
            xaxis: {...defaultLayout.xaxis, title: 'Minimum Coverage Depth', range: [0, 100]},
            yaxis: {...defaultLayout.yaxis, title: '% of Genome', range: [0, 105]},
        }, plotConfig);
    })();
    {% endif %}

    {% if has_gc_bias %}
    // GC Bias
    (function() {
        const traces = DATA.filter(s => s.gc_bias_coverage_by_gc).map(s => {
            // Normalize to mean coverage
            const meanCov = s.genome_mean_coverage || 1;
            const normalized = s.gc_bias_coverage_by_gc.map(v => v / meanCov);
            const x = normalized.map((_, i) => i);
            return {
                x: x,
                y: normalized,
                type: 'scatter',
                mode: 'lines',
                name: s.sample_name,
                line: {color: s.color, width: 2},
            };
        });
        // Add ideal line at y=1
        traces.push({
            x: [0, 100],
            y: [1, 1],
            type: 'scatter',
            mode: 'lines',
            name: 'Ideal',
            line: {color: '#999', width: 1, dash: 'dash'},
            showlegend: true,
        });
        Plotly.newPlot('plot-gc-bias', traces, {
            ...defaultLayout,
            xaxis: {...defaultLayout.xaxis, title: 'GC Content (%)'},
            yaxis: {...defaultLayout.yaxis, title: 'Normalized Coverage'},
        }, plotConfig);
    })();
    {% endif %}

    // Mapping Quality
    (function() {
        const traces = DATA.map(s => ({
            x: s.mapq_hist_keys,
            y: s.mapq_hist_values,
            type: 'bar',
            name: s.sample_name,
            marker: {color: s.color},
            opacity: 0.7,
        }));
        Plotly.newPlot('plot-mapq', traces, {
            ...defaultLayout,
            barmode: 'overlay',
            xaxis: {...defaultLayout.xaxis, title: 'Mapping Quality (MAPQ)'},
            yaxis: {...defaultLayout.yaxis, title: 'Count'},
        }, plotConfig);
    })();

    // Indels
    (function() {
        // Insertion size distribution
        const traces = [];
        DATA.forEach(s => {
            if (s.insertion_size_keys && s.insertion_size_keys.length > 0) {
                traces.push({
                    x: s.insertion_size_keys,
                    y: s.insertion_size_values,
                    type: 'bar',
                    name: s.sample_name + ' (ins)',
                    marker: {color: s.color},
                    opacity: 0.6,
                });
            }
            if (s.deletion_size_keys && s.deletion_size_keys.length > 0) {
                traces.push({
                    x: s.deletion_size_keys.map(v => -v),  // negative x for deletions
                    y: s.deletion_size_values,
                    type: 'bar',
                    name: s.sample_name + ' (del)',
                    marker: {color: s.color, opacity: 0.4},
                });
            }
        });
        Plotly.newPlot('plot-indels', traces, {
            ...defaultLayout,
            barmode: 'overlay',
            xaxis: {...defaultLayout.xaxis, title: 'Size (bp) — negative = deletions, positive = insertions'},
            yaxis: {...defaultLayout.yaxis, title: 'Count'},
        }, plotConfig);
    })();

    // Sidebar active state tracking
    (function() {
        const sections = document.querySelectorAll('.section');
        const links = document.querySelectorAll('.sidebar a');
        window.addEventListener('scroll', () => {
            let current = '';
            sections.forEach(s => {
                if (window.scrollY >= s.offsetTop - 100) current = s.id;
            });
            links.forEach(a => {
                a.classList.toggle('active', a.getAttribute('href') === '#' + current);
            });
        });
    })();
    </script>
</body>
</html>
```

**Step 3: Update python/hifi_qc/cli.py to call report generator**

```python
import click


@click.command()
@click.argument("bam_files", nargs=-1, required=True, type=click.Path(exists=True))
@click.option("-o", "--output", default="hifi_qc_report.html", help="Output HTML report path")
@click.option("-r", "--reference", default=None, type=click.Path(exists=True), help="Reference FASTA for GC bias")
@click.option("-t", "--threads", default=None, type=int, help="Number of threads")
@click.option("--sample-fraction", default=None, type=float, help="Fraction of reads to sample (0.0-1.0)")
@click.option("--seed", default=None, type=int, help="RNG seed for reproducible sampling")
@click.option("--sample-names", default=None, type=str, help="Comma-separated sample names")
def main(bam_files, output, reference, threads, sample_fraction, seed, sample_names):
    """Quality control for PacBio HiFi sequencing data."""
    from hifi_qc._core import process_bam_files
    from hifi_qc.report import generate_report

    names = sample_names.split(",") if sample_names else None

    click.echo(f"Processing {len(bam_files)} file(s)...")
    results = process_bam_files(
        list(bam_files),
        reference=reference,
        threads=threads,
        sample_fraction=sample_fraction,
        seed=seed,
        sample_names=names,
    )

    click.echo(f"Generating report for {len(results)} sample(s)...")
    generate_report(results, output)
    click.echo(f"Report written to: {output}")
```

**Step 4: Build and verify**

Run: `cd /Users/arq5x/src/quinlan/hifi-qc && maturin develop`
Expected: Builds.

**Step 5: Commit**

```bash
git add python/hifi_qc/report.py python/hifi_qc/templates/report.html.j2 python/hifi_qc/cli.py
git commit -m "feat: add HTML report template with Plotly.js charts"
```

---

## Task 8: Create test BAM file and end-to-end test

**Files:**
- Create: `tests/create_test_bam.py`
- Create: `tests/test_e2e.py`

We need a small synthetic BAM file with HiFi-like reads (np/rq tags, realistic CIGAR ops) for testing.

**Step 1: Create tests/create_test_bam.py**

This script uses pysam (install separately for test data generation only) to create a small test BAM file with:
- A tiny reference (2 chromosomes, 10kb each)
- ~100 reads with realistic HiFi properties
- np and rq tags set
- Mix of M, I, D, =, X CIGAR ops

```python
"""Generate a small synthetic HiFi BAM file for testing."""

import pysam
import random
import os

random.seed(42)

CHROMS = [("chr1", 10000), ("chr2", 8000)]
N_READS = 100
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


def random_sequence(length):
    return "".join(random.choice("ACGT") for _ in range(length))


def create_test_bam():
    bam_path = os.path.join(OUTPUT_DIR, "data", "test.bam")
    ref_path = os.path.join(OUTPUT_DIR, "data", "test_ref.fa")
    os.makedirs(os.path.join(OUTPUT_DIR, "data"), exist_ok=True)

    # Create reference FASTA
    with open(ref_path, "w") as f:
        for name, length in CHROMS:
            seq = random_sequence(length)
            f.write(f">{name}\n{seq}\n")

    # Index it
    pysam.faidx(ref_path)

    # Create BAM header
    header = pysam.AlignmentHeader.from_references(
        [name for name, _ in CHROMS],
        [length for _, length in CHROMS],
    )

    # Write BAM
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for i in range(N_READS):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i:04d}"

            # Random chromosome
            chrom_idx = random.randint(0, len(CHROMS) - 1)
            chrom_name, chrom_len = CHROMS[chrom_idx]
            a.reference_id = chrom_idx

            # Random read length (HiFi-like: 5-20kb)
            read_len = random.randint(5000, 20000)

            # Random position (ensure it fits)
            max_pos = chrom_len - read_len - 100
            if max_pos < 0:
                max_pos = 0
            a.reference_start = random.randint(0, max_pos)

            # Generate sequence and quality
            a.query_sequence = random_sequence(read_len)
            a.query_qualities = pysam.qualitystring_to_array(
                "".join(chr(random.randint(30, 40) + 33) for _ in range(read_len))
            )

            # Simple CIGAR: mostly matches with occasional small indels
            cigar = []
            remaining = read_len
            while remaining > 0:
                op_len = min(random.randint(100, 2000), remaining)
                cigar.append((0, op_len))  # M
                remaining -= op_len
                if remaining > 5 and random.random() < 0.1:
                    ins_len = random.randint(1, 3)
                    cigar.append((1, ins_len))  # I
                    remaining -= ins_len
                if remaining > 5 and random.random() < 0.1:
                    cigar.append((2, random.randint(1, 3)))  # D
            a.cigar = cigar

            a.mapping_quality = 60
            a.flag = 0  # mapped, primary

            # HiFi tags
            a.set_tag("np", random.randint(3, 30), "i")  # number of passes
            a.set_tag("rq", round(random.uniform(0.99, 0.9999), 4), "f")  # read quality

            bam.write(a)

    # Sort and index
    sorted_path = bam_path.replace(".bam", ".sorted.bam")
    pysam.sort("-o", sorted_path, bam_path)
    os.replace(sorted_path, bam_path)
    pysam.index(bam_path)

    print(f"Created test BAM: {bam_path}")
    print(f"Created test reference: {ref_path}")


if __name__ == "__main__":
    create_test_bam()
```

**Step 2: Create tests/test_e2e.py**

```python
"""End-to-end test for hifi-qc."""

import os
import subprocess
import sys

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(TEST_DIR, "data")
BAM_FILE = os.path.join(DATA_DIR, "test.bam")
REF_FILE = os.path.join(DATA_DIR, "test_ref.fa")
OUTPUT_FILE = os.path.join(DATA_DIR, "test_report.html")


def test_basic_run():
    """Test basic hifi-qc run with a BAM file."""
    result = subprocess.run(
        ["hifi-qc", BAM_FILE, "-o", OUTPUT_FILE],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"hifi-qc failed: {result.stderr}"
    assert os.path.exists(OUTPUT_FILE), "Report file not created"
    with open(OUTPUT_FILE) as f:
        html = f.read()
    assert "hifi-qc Report" in html
    assert "General Statistics" in html
    assert "Read Length" in html
    print("PASS: basic run")


def test_with_reference():
    """Test hifi-qc with reference FASTA for GC bias."""
    output = OUTPUT_FILE.replace(".html", "_ref.html")
    result = subprocess.run(
        ["hifi-qc", BAM_FILE, "-o", output, "-r", REF_FILE],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"hifi-qc failed: {result.stderr}"
    with open(output) as f:
        html = f.read()
    assert "GC Bias" in html
    print("PASS: with reference")


def test_sampling_mode():
    """Test hifi-qc with sampling mode."""
    output = OUTPUT_FILE.replace(".html", "_sampled.html")
    result = subprocess.run(
        ["hifi-qc", BAM_FILE, "-o", output, "--sample-fraction", "0.5", "--seed", "42"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"hifi-qc failed: {result.stderr}"
    with open(output) as f:
        html = f.read()
    # Coverage sections should be absent in sampling mode
    assert "Sampling mode" in html or "sampling" in html.lower()
    print("PASS: sampling mode")


if __name__ == "__main__":
    test_basic_run()
    test_with_reference()
    test_sampling_mode()
    print("\nAll tests passed!")
```

**Step 3: Generate test data and run tests**

Run: `pip install pysam && cd /Users/arq5x/src/quinlan/hifi-qc && python tests/create_test_bam.py`
Expected: Creates `tests/data/test.bam` and `tests/data/test_ref.fa`.

Run: `cd /Users/arq5x/src/quinlan/hifi-qc && python tests/test_e2e.py`
Expected: All tests pass.

**Step 4: Commit**

```bash
git add tests/
git commit -m "test: add synthetic BAM generation and end-to-end tests"
```

---

## Task 9: Iteration and polish

**Files:**
- Potentially any Rust or Python file

This is the debugging/polish phase. After the end-to-end test runs, review the output HTML report and fix issues:

**Step 1: Open the report in a browser and check each section**

Run: `open /Users/arq5x/src/quinlan/hifi-qc/tests/data/test_report.html`

Check:
- Summary table renders with correct values
- All histograms have data and look reasonable
- Scatter plot renders
- Coverage charts (if coverage test run) show data
- Sidebar navigation works
- Charts are interactive (zoom, hover)

**Step 2: Fix any Rust compilation issues or runtime errors**

Common issues to watch for:
- noodles API mismatches (method names, type conversions) — fix based on compiler errors
- PyO3 type conversion issues — ensure all HashMap keys/values are convertible
- Integer overflow in depth arrays — `saturating_add` handles this
- Empty data edge cases — check for division by zero in stats

**Step 3: Fix any HTML/JS rendering issues**

Common issues:
- Plotly.js CDN not loading — verify script tag URL
- JSON serialization of large arrays — check for NaN/Infinity values
- Layout overflow — verify responsive CSS

**Step 4: Commit fixes**

```bash
git add -A
git commit -m "fix: address issues found during end-to-end testing"
```

---

## Task 10: Final integration test with real data (manual)

**Step 1: Test with a real HiFi BAM file (if available)**

Run: `hifi-qc /path/to/real/hifi.bam -o real_report.html -r /path/to/ref.fa -t 8`

Check:
- Completes without errors
- Report looks reasonable for real data
- Coverage values match expectations (~30x for standard WGS)
- Performance is acceptable (should process ~30x human WGS in < 30 minutes)

**Step 2: Test multi-sample mode**

Run: `hifi-qc sample1.bam sample2.bam -o multi_report.html`

Check:
- Both samples appear in summary table
- Charts show overlaid data with different colors

**Step 3: Commit any final fixes**

```bash
git add -A
git commit -m "fix: address issues from real data testing"
```
