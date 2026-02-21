pub mod bam_reader;

use pyo3::prelude::*;

/// Placeholder: process BAM/CRAM files and return QC metrics.
#[pyfunction]
#[pyo3(signature = (bam_paths, reference=None, threads=None, sample_fraction=None, seed=None))]
fn process_bam_files(
    bam_paths: Vec<String>,
    reference: Option<String>,
    threads: Option<usize>,
    sample_fraction: Option<f64>,
    seed: Option<u64>,
) -> PyResult<Vec<PyObject>> {
    // Suppress unused variable warnings for now
    let _ = (bam_paths, reference, threads, sample_fraction, seed);
    Ok(vec![])
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(process_bam_files, m)?)?;
    Ok(())
}
