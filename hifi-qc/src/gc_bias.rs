//! GC bias computation module.
//!
//! Reads a reference FASTA in non-overlapping windows and stratifies observed
//! per-base coverage by GC content to quantify GC bias.

use std::io;
use std::path::Path;

use noodles::fasta;

/// Window size (in bases) used to tile the reference genome.
const WINDOW_SIZE: usize = 100;

/// Number of GC bins (0% through 100%, inclusive).
const NUM_GC_BINS: usize = 101;

/// Result of GC bias computation.
pub struct GcBiasResult {
    /// For each GC bin (0..=100): number of reference windows in that bin.
    pub reference_gc_distribution: Vec<u64>,
    /// For each GC bin (0..=100): mean coverage of windows in that bin.
    pub coverage_by_gc: Vec<f64>,
}

/// Compute GC bias by tiling the reference genome in [`WINDOW_SIZE`]-bp windows
/// and stratifying the observed depth-of-coverage by GC content.
///
/// # Arguments
///
/// * `reference_path` - Path to an indexed FASTA file (`.fai` must exist).
/// * `depth_arrays`   - One depth array per reference sequence, in BAM header order.
///                      Each inner `Vec<u32>` has one entry per base.
/// * `reference_names` - Chromosome / contig names in BAM header order.
///
/// # Errors
///
/// Returns an `io::Error` if the FASTA file cannot be opened or a chromosome
/// cannot be queried.
pub fn compute_gc_bias<P: AsRef<Path>>(
    reference_path: P,
    depth_arrays: &[Vec<u32>],
    reference_names: &[String],
) -> io::Result<GcBiasResult> {
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference_path)?;

    // Accumulators: sum of mean-coverage values and window counts per GC bin.
    let mut gc_coverage_sums = vec![0.0_f64; NUM_GC_BINS];
    let mut gc_window_counts = vec![0_u64; NUM_GC_BINS];

    for (chrom_idx, chrom_name) in reference_names.iter().enumerate() {
        // If we have no depth array for this chromosome, skip it.
        let depth = match depth_arrays.get(chrom_idx) {
            Some(d) => d,
            None => continue,
        };

        // Query the full chromosome sequence from the indexed FASTA.
        let region: noodles::core::Region = chrom_name
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("{}", e)))?;

        let record = reader.query(&region)?;
        let seq: Vec<u8> = record.sequence().as_ref().iter().copied().collect();

        let chrom_len = seq.len();

        // Tile the chromosome in non-overlapping windows.
        let num_full_windows = chrom_len / WINDOW_SIZE;
        for win_idx in 0..num_full_windows {
            let win_start = win_idx * WINDOW_SIZE;
            let win_end = win_start + WINDOW_SIZE;
            let window_seq = &seq[win_start..win_end];

            // Count GC and N bases in this window.
            let mut gc_count: usize = 0;
            let mut n_count: usize = 0;
            for &base in window_seq {
                match base {
                    b'G' | b'g' | b'C' | b'c' => gc_count += 1,
                    b'N' | b'n' => n_count += 1,
                    _ => {}
                }
            }

            // Skip windows where N bases exceed 50% of the window.
            if n_count > WINDOW_SIZE / 2 {
                continue;
            }

            // GC% = gc_count / (window_size - n_count) * 100, rounded to nearest int.
            let effective_len = WINDOW_SIZE - n_count;
            let gc_pct = if effective_len > 0 {
                (gc_count as f64 / effective_len as f64 * 100.0).round() as usize
            } else {
                continue; // All N's -- skip
            };
            let gc_bin = gc_pct.min(100);

            // Mean coverage across the window from the depth array.
            let depth_slice_end = win_end.min(depth.len());
            let depth_slice_start = win_start.min(depth.len());
            if depth_slice_start >= depth_slice_end {
                continue;
            }
            let mean_cov: f64 = depth[depth_slice_start..depth_slice_end]
                .iter()
                .map(|&d| d as f64)
                .sum::<f64>()
                / WINDOW_SIZE as f64;

            gc_coverage_sums[gc_bin] += mean_cov;
            gc_window_counts[gc_bin] += 1;
        }
    }

    // Compute final per-bin mean coverage.
    let coverage_by_gc: Vec<f64> = gc_coverage_sums
        .iter()
        .zip(gc_window_counts.iter())
        .map(|(&sum, &count)| {
            if count > 0 {
                sum / count as f64
            } else {
                0.0
            }
        })
        .collect();

    Ok(GcBiasResult {
        reference_gc_distribution: gc_window_counts,
        coverage_by_gc,
    })
}
