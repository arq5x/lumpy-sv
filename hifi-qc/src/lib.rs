pub mod bam_reader;
pub mod coverage;
pub mod gc_bias;
pub mod metrics;

use std::path::Path;

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// Helper: collect a HashMap into two sorted Vecs (keys, values), sorted by key.
fn sorted_hist_u32(map: &std::collections::HashMap<u32, u64>) -> (Vec<u32>, Vec<u64>) {
    let mut entries: Vec<(u32, u64)> = map.iter().map(|(&k, &v)| (k, v)).collect();
    entries.sort_by_key(|&(k, _)| k);
    entries.into_iter().unzip()
}

fn sorted_hist_u8(map: &std::collections::HashMap<u8, u64>) -> (Vec<u8>, Vec<u64>) {
    let mut entries: Vec<(u8, u64)> = map.iter().map(|(&k, &v)| (k, v)).collect();
    entries.sort_by_key(|&(k, _)| k);
    entries.into_iter().unzip()
}

/// Build a Python dict from the accumulated metrics and optional coverage/gc_bias results.
fn build_sample_dict(
    py: Python<'_>,
    m: &metrics::SampleMetrics,
    header: &bam_reader::FileHeader,
    coverage_opt: Option<&coverage::CoverageSummary>,
    gc_bias_opt: Option<&gc_bias::GcBiasResult>,
) -> PyResult<PyObject> {
    let dict = PyDict::new(py);

    // Basic counts
    dict.set_item("sample_name", &m.sample_name)?;
    dict.set_item("total_reads", m.total_reads)?;
    dict.set_item("mapped_reads", m.mapped_reads)?;
    dict.set_item("unmapped_reads", m.unmapped_reads)?;
    dict.set_item("duplicate_reads", m.duplicate_reads)?;
    dict.set_item("total_bases", m.total_bases)?;

    // Summary statistics
    dict.set_item("mean_read_length", m.mean_read_length())?;
    dict.set_item("median_read_length", m.median_read_length())?;
    dict.set_item("n50_read_length", m.compute_n50())?;
    dict.set_item("mean_qscore", m.mean_qscore())?;

    // Sampling flag
    dict.set_item("is_sampled", m.is_sampled)?;

    // Read length histogram (sorted by key)
    let (len_keys, len_vals) = sorted_hist_u32(&m.length_histogram);
    dict.set_item("length_hist_keys", len_keys)?;
    dict.set_item("length_hist_values", len_vals)?;

    // Q-score histogram (sorted by key)
    let (qs_keys, qs_vals) = sorted_hist_u8(&m.qscore_histogram);
    dict.set_item("qscore_hist_keys", qs_keys)?;
    dict.set_item("qscore_hist_values", qs_vals)?;

    // Length-quality pairs
    let (lq_lengths, lq_scores): (Vec<u32>, Vec<f32>) =
        m.length_quality_pairs.iter().copied().unzip();
    dict.set_item("length_quality_lengths", lq_lengths)?;
    dict.set_item("length_quality_scores", lq_scores)?;

    // HiFi passes histogram (sorted by key)
    let (passes_keys, passes_vals) = sorted_hist_u32(&m.passes_histogram);
    dict.set_item("passes_hist_keys", passes_keys)?;
    dict.set_item("passes_hist_values", passes_vals)?;

    // Read quality histogram (sorted by key, keys converted to f64 / 10000)
    let (rq_keys_raw, rq_vals) = sorted_hist_u32(&m.rq_histogram);
    let rq_keys_f64: Vec<f64> = rq_keys_raw.iter().map(|&k| k as f64 / 10000.0).collect();
    dict.set_item("rq_hist_keys", rq_keys_f64)?;
    dict.set_item("rq_hist_values", rq_vals)?;

    // Mapping quality histogram (sorted by key)
    let (mapq_keys, mapq_vals) = sorted_hist_u8(&m.mapq_histogram);
    dict.set_item("mapq_hist_keys", mapq_keys)?;
    dict.set_item("mapq_hist_values", mapq_vals)?;

    // GC content histogram (sorted by key)
    let (gc_keys, gc_vals) = sorted_hist_u8(&m.gc_content_histogram);
    dict.set_item("gc_content_keys", gc_keys)?;
    dict.set_item("gc_content_values", gc_vals)?;

    // Mismatch and indel statistics
    dict.set_item("total_mismatches", m.total_mismatches)?;
    dict.set_item("total_insertions", m.total_insertions)?;
    dict.set_item("total_insertion_bases", m.total_insertion_bases)?;
    dict.set_item("total_deletions", m.total_deletions)?;
    dict.set_item("total_deletion_bases", m.total_deletion_bases)?;

    // Insertion size histogram (sorted by key)
    let (ins_keys, ins_vals) = sorted_hist_u32(&m.insertion_size_histogram);
    dict.set_item("insertion_size_keys", ins_keys)?;
    dict.set_item("insertion_size_values", ins_vals)?;

    // Deletion size histogram (sorted by key)
    let (del_keys, del_vals) = sorted_hist_u32(&m.deletion_size_histogram);
    dict.set_item("deletion_size_keys", del_keys)?;
    dict.set_item("deletion_size_values", del_vals)?;

    // Coverage (only if not sampling)
    if let Some(cov) = coverage_opt {
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

    // GC bias (only if not sampling and reference provided)
    if let Some(gc) = gc_bias_opt {
        dict.set_item("gc_bias_ref_distribution", &gc.reference_gc_distribution)?;
        dict.set_item("gc_bias_coverage_by_gc", &gc.coverage_by_gc)?;
    }

    // Also include header info for reference
    let _ = header; // header data already incorporated via coverage

    Ok(dict.into())
}

/// Process BAM/CRAM files and return QC metrics as a list of Python dicts.
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
    // Set up rayon thread pool if threads specified
    if let Some(n) = threads {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global();
        // Ignore error if pool already initialized (e.g., called twice)
    }

    let is_sampling = sample_fraction.is_some();
    let fraction = sample_fraction.unwrap_or(1.0);
    let base_seed = seed.unwrap_or(42);
    let ref_str = reference.as_deref();

    let mut all_results: Vec<PyObject> = Vec::with_capacity(bam_paths.len());

    for (file_idx, bam_path) in bam_paths.iter().enumerate() {
        // (a) Read header
        let header = bam_reader::read_header(bam_path, ref_str)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!(
                "Failed to read header from {}: {}", bam_path, e
            )))?;

        // (b) Determine sample name
        let sample_name = if let Some(ref names) = sample_names {
            if file_idx < names.len() {
                names[file_idx].clone()
            } else {
                // Fallback to header or filename
                header.sample_name.clone().unwrap_or_else(|| {
                    Path::new(bam_path)
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("unknown")
                        .to_string()
                })
            }
        } else {
            // Try header RG sample name, then filename stem
            header.sample_name.clone().unwrap_or_else(|| {
                Path::new(bam_path)
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .to_string()
            })
        };

        // (c) Create metrics accumulator
        let mut m = metrics::SampleMetrics::new(
            sample_name,
            &header.reference_lengths,
            is_sampling,
        );
        m.is_sampled = is_sampling;

        // (d) Set up optional RNG for sampling
        let mut rng: Option<StdRng> = if is_sampling {
            Some(StdRng::seed_from_u64(base_seed + file_idx as u64))
        } else {
            None
        };

        // (e) Process alignment file with closure
        bam_reader::process_alignment_file(bam_path, ref_str, |read| {
            // If sampling, skip reads randomly
            if let Some(ref mut r) = rng {
                if r.random::<f64>() >= fraction {
                    return;
                }
            }
            m.process_read(&read);
        })
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!(
            "Failed to process {}: {}", bam_path, e
        )))?;

        // (f) Summarize coverage (skip if sampling)
        let coverage_opt = if !is_sampling && !m.depth_arrays.is_empty() {
            Some(coverage::summarize_coverage(
                &m.depth_arrays,
                &header.reference_names,
                1000,
            ))
        } else {
            None
        };

        // (g) Compute GC bias (skip if sampling or no reference)
        let gc_bias_opt = if !is_sampling {
            if let Some(ref ref_path) = reference {
                match gc_bias::compute_gc_bias(ref_path, &m.depth_arrays, &header.reference_names) {
                    Ok(result) => Some(result),
                    Err(_) => None, // Silently skip if GC bias computation fails
                }
            } else {
                None
            }
        } else {
            None
        };

        // (h) Build Python dict with all results
        let dict = build_sample_dict(
            py,
            &m,
            &header,
            coverage_opt.as_ref(),
            gc_bias_opt.as_ref(),
        )?;
        all_results.push(dict);
    }

    // Return as a Python list
    let list = PyList::new(py, &all_results)?;
    Ok(list.into())
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(process_bam_files, m)?)?;
    Ok(())
}
