//! Coverage summarization module.
//!
//! Takes populated depth arrays from `SampleMetrics` and computes coverage
//! statistics in parallel across chromosomes using rayon.

use rayon::prelude::*;

/// Per-chromosome coverage statistics.
pub struct ChromCoverage {
    pub name: String,
    pub length: usize,
    pub mean: f64,
    pub median: f64,
    pub stdev: f64,
}

/// Genome-wide coverage summary computed from per-base depth arrays.
pub struct CoverageSummary {
    pub per_chrom: Vec<ChromCoverage>,
    pub genome_mean: f64,
    pub genome_median: f64,
    pub genome_stdev: f64,
    /// Index = depth, value = number of bases at that depth.
    pub coverage_histogram: Vec<u64>,
    /// Fraction of genome covered at >= each depth threshold (index 0..max_depth).
    pub cumulative_coverage: Vec<f64>,
}

/// Compute the median from a histogram of counts.
///
/// Walks bins in order, accumulating counts until reaching the position
/// that corresponds to the 50th percentile.
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

/// Build a depth histogram for a single chromosome's depth array.
///
/// Depths above `max_depth - 1` are capped into the last bin.
fn build_chrom_histogram(depths: &[u32], max_depth: usize) -> Vec<u64> {
    let mut hist = vec![0u64; max_depth];
    for &d in depths {
        let bin = (d as usize).min(max_depth - 1);
        hist[bin] += 1;
    }
    hist
}

/// Compute per-chromosome statistics from its depth array.
///
/// Returns `(ChromCoverage, histogram)` so that the histogram can be
/// merged into the genome-wide histogram without a second pass.
fn compute_chrom_stats(name: &str, depths: &[u32], max_depth: usize) -> (ChromCoverage, Vec<u64>) {
    let length = depths.len();
    if length == 0 {
        let cov = ChromCoverage {
            name: name.to_string(),
            length: 0,
            mean: 0.0,
            median: 0.0,
            stdev: 0.0,
        };
        return (cov, vec![0u64; max_depth]);
    }

    let hist = build_chrom_histogram(depths, max_depth);
    let mean = {
        let sum: u64 = hist
            .iter()
            .enumerate()
            .map(|(depth, &count)| depth as u64 * count)
            .sum();
        sum as f64 / length as f64
    };
    let median = histogram_median(&hist, length);
    let stdev = {
        let variance: f64 = hist
            .iter()
            .enumerate()
            .map(|(depth, &count)| {
                let diff = depth as f64 - mean;
                diff * diff * count as f64
            })
            .sum::<f64>()
            / length as f64;
        variance.sqrt()
    };

    let cov = ChromCoverage {
        name: name.to_string(),
        length,
        mean,
        median,
        stdev,
    };
    (cov, hist)
}

/// Summarize coverage across all chromosomes in parallel.
///
/// # Arguments
///
/// * `depth_arrays` - One depth array per reference sequence (from `SampleMetrics`).
/// * `reference_names` - Names corresponding to each depth array.
/// * `max_depth` - Cap for the coverage histogram. Depths at or above this
///   value are placed in the last bin (index `max_depth - 1`).
pub fn summarize_coverage(
    depth_arrays: &[Vec<u32>],
    reference_names: &[String],
    max_depth: usize,
) -> CoverageSummary {
    let max_depth = max_depth.max(1); // ensure at least 1 bin

    // Compute per-chromosome stats in parallel using rayon.
    let chrom_results: Vec<(ChromCoverage, Vec<u64>)> = depth_arrays
        .par_iter()
        .enumerate()
        .map(|(i, depths)| {
            let name = reference_names
                .get(i)
                .map(|s| s.as_str())
                .unwrap_or("unknown");
            compute_chrom_stats(name, depths, max_depth)
        })
        .collect();

    // Separate ChromCoverage structs from their histograms.
    let mut per_chrom = Vec::with_capacity(chrom_results.len());
    let mut genome_histogram = vec![0u64; max_depth];

    for (cov, hist) in chrom_results {
        for (bin, &count) in hist.iter().enumerate() {
            genome_histogram[bin] += count;
        }
        per_chrom.push(cov);
    }

    // Genome-wide total bases.
    let total_bases: u64 = genome_histogram.iter().sum();

    // Genome-wide mean from histogram.
    let genome_mean = if total_bases == 0 {
        0.0
    } else {
        let weighted_sum: u64 = genome_histogram
            .iter()
            .enumerate()
            .map(|(depth, &count)| depth as u64 * count)
            .sum();
        weighted_sum as f64 / total_bases as f64
    };

    // Genome-wide median from histogram.
    let genome_median = histogram_median(&genome_histogram, total_bases as usize);

    // Genome-wide stdev from histogram.
    let genome_stdev = if total_bases == 0 {
        0.0
    } else {
        let variance: f64 = genome_histogram
            .iter()
            .enumerate()
            .map(|(depth, &count)| {
                let diff = depth as f64 - genome_mean;
                diff * diff * count as f64
            })
            .sum::<f64>()
            / total_bases as f64;
        variance.sqrt()
    };

    // Cumulative coverage: fraction of genome at >= each depth threshold.
    let cumulative_coverage = if total_bases == 0 {
        vec![0.0; max_depth]
    } else {
        // Build a suffix sum of the histogram.
        let mut suffix_sum = vec![0u64; max_depth];
        suffix_sum[max_depth - 1] = genome_histogram[max_depth - 1];
        for i in (0..max_depth - 1).rev() {
            suffix_sum[i] = suffix_sum[i + 1] + genome_histogram[i];
        }
        suffix_sum
            .iter()
            .map(|&count| count as f64 / total_bases as f64)
            .collect()
    };

    CoverageSummary {
        per_chrom,
        genome_mean,
        genome_median,
        genome_stdev,
        coverage_histogram: genome_histogram,
        cumulative_coverage,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_histogram_median_empty() {
        assert_eq!(histogram_median(&[], 0), 0.0);
        assert_eq!(histogram_median(&[0, 0, 0], 0), 0.0);
    }

    #[test]
    fn test_histogram_median_simple() {
        // 10 bases at depth 0, 20 at depth 1, 10 at depth 2
        // total = 40, half = 20
        // cumulative: bin0 = 10 (<=20), bin1 = 30 (>20) -> median = 1.0
        let hist = vec![10u64, 20, 10];
        assert_eq!(histogram_median(&hist, 40), 1.0);
    }

    #[test]
    fn test_histogram_median_all_zero_depth() {
        // All 100 bases at depth 0
        let hist = vec![100u64, 0, 0];
        // total=100, half=50. cumulative after bin0 = 100 > 50 -> median = 0
        assert_eq!(histogram_median(&hist, 100), 0.0);
    }

    #[test]
    fn test_histogram_median_all_same_depth() {
        // All 50 bases at depth 5
        let mut hist = vec![0u64; 10];
        hist[5] = 50;
        assert_eq!(histogram_median(&hist, 50), 5.0);
    }

    #[test]
    fn test_empty_depth_arrays() {
        let result = summarize_coverage(&[], &[], 1000);
        assert!(result.per_chrom.is_empty());
        assert_eq!(result.genome_mean, 0.0);
        assert_eq!(result.genome_median, 0.0);
        assert_eq!(result.genome_stdev, 0.0);
        assert_eq!(result.coverage_histogram.len(), 1000);
        assert!(result.coverage_histogram.iter().all(|&c| c == 0));
        assert_eq!(result.cumulative_coverage.len(), 1000);
        assert!(result.cumulative_coverage.iter().all(|&f| f == 0.0));
    }

    #[test]
    fn test_single_chrom_uniform_depth() {
        // 100 bases all at depth 10
        let depths = vec![10u32; 100];
        let names = vec!["chr1".to_string()];
        let result = summarize_coverage(&[depths], &names, 1000);

        assert_eq!(result.per_chrom.len(), 1);
        let chrom = &result.per_chrom[0];
        assert_eq!(chrom.name, "chr1");
        assert_eq!(chrom.length, 100);
        assert!((chrom.mean - 10.0).abs() < 1e-10);
        assert!((chrom.median - 10.0).abs() < 1e-10);
        assert!(chrom.stdev.abs() < 1e-10); // no variance

        assert!((result.genome_mean - 10.0).abs() < 1e-10);
        assert!((result.genome_median - 10.0).abs() < 1e-10);
        assert!(result.genome_stdev.abs() < 1e-10);

        // Histogram: all 100 bases in bin 10
        assert_eq!(result.coverage_histogram[10], 100);
        assert_eq!(result.coverage_histogram[0], 0);
        assert_eq!(result.coverage_histogram[11], 0);

        // Cumulative: at depth 0..=10, should be 1.0; at depth 11, should be 0.0
        assert!((result.cumulative_coverage[0] - 1.0).abs() < 1e-10);
        assert!((result.cumulative_coverage[10] - 1.0).abs() < 1e-10);
        assert!((result.cumulative_coverage[11] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_single_chrom_varying_depth() {
        // 50 bases at depth 0, 50 bases at depth 20
        let mut depths = vec![0u32; 50];
        depths.extend(vec![20u32; 50]);
        let names = vec!["chrX".to_string()];
        let result = summarize_coverage(&[depths], &names, 1000);

        let chrom = &result.per_chrom[0];
        assert_eq!(chrom.length, 100);
        // mean = (0*50 + 20*50) / 100 = 10.0
        assert!((chrom.mean - 10.0).abs() < 1e-10);
        // median: total=100, half=50. cumulative: bin0=50 (<=50), bin20=100 (>50) -> 20
        // Wait: half=50, cumulative after bin0 = 50, which is NOT > 50. So we keep going.
        // bin1..19 are 0, cumulative stays 50. bin20 = 100 > 50 -> median = 20
        assert!((chrom.median - 20.0).abs() < 1e-10);

        // stdev = sqrt((0-10)^2 * 50 + (20-10)^2 * 50) / 100) = sqrt(100*50+100*50)/100 = sqrt(100) = 10
        assert!((chrom.stdev - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_multiple_chroms_parallel() {
        // chr1: 100 bases at depth 5
        // chr2: 200 bases at depth 10
        let depths1 = vec![5u32; 100];
        let depths2 = vec![10u32; 200];
        let names = vec!["chr1".to_string(), "chr2".to_string()];
        let result = summarize_coverage(&[depths1, depths2], &names, 1000);

        assert_eq!(result.per_chrom.len(), 2);

        // Find chr1 and chr2 (order should be preserved by rayon indexed iteration)
        let chr1 = &result.per_chrom[0];
        let chr2 = &result.per_chrom[1];

        assert_eq!(chr1.name, "chr1");
        assert_eq!(chr1.length, 100);
        assert!((chr1.mean - 5.0).abs() < 1e-10);

        assert_eq!(chr2.name, "chr2");
        assert_eq!(chr2.length, 200);
        assert!((chr2.mean - 10.0).abs() < 1e-10);

        // Genome-wide: total = 300 bases
        // weighted sum = 5*100 + 10*200 = 500 + 2000 = 2500
        // mean = 2500/300 = 8.333...
        let expected_mean = 2500.0 / 300.0;
        assert!((result.genome_mean - expected_mean).abs() < 1e-10);

        // Histogram: bin 5 = 100, bin 10 = 200
        assert_eq!(result.coverage_histogram[5], 100);
        assert_eq!(result.coverage_histogram[10], 200);

        // Genome median: total=300, half=150
        // cumulative: bins 0..4 = 0, bin5 = 100 (<=150), bins 6..9 = 100, bin10 = 300 > 150
        // median = 10
        assert!((result.genome_median - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_max_depth_capping() {
        // Depths exceeding max_depth should go into the last bin
        let depths = vec![999u32; 50]; // all at depth 999
        let names = vec!["chr1".to_string()];
        let max_depth = 100; // cap at 100

        let result = summarize_coverage(&[depths], &names, max_depth);

        // Depth 999 should be capped into bin 99 (max_depth - 1)
        assert_eq!(result.coverage_histogram[99], 50);
        assert_eq!(result.coverage_histogram.len(), 100);

        // Mean from histogram perspective: all in bin 99, so mean = 99
        // (the actual depths are 999, but histogram-based stats use the capped value)
        assert!((result.genome_mean - 99.0).abs() < 1e-10);
    }

    #[test]
    fn test_cumulative_coverage() {
        // 50 bases at depth 0, 30 bases at depth 5, 20 bases at depth 10
        let mut depths = vec![0u32; 50];
        depths.extend(vec![5u32; 30]);
        depths.extend(vec![10u32; 20]);
        let names = vec!["chr1".to_string()];
        let result = summarize_coverage(&[depths], &names, 20);

        // Total = 100 bases
        // cumulative_coverage[0] = fraction at >= 0 = 100/100 = 1.0
        assert!((result.cumulative_coverage[0] - 1.0).abs() < 1e-10);
        // cumulative_coverage[1] = fraction at >= 1 = (30 + 20) / 100 = 0.5
        assert!((result.cumulative_coverage[1] - 0.5).abs() < 1e-10);
        // cumulative_coverage[5] = fraction at >= 5 = (30 + 20) / 100 = 0.5
        assert!((result.cumulative_coverage[5] - 0.5).abs() < 1e-10);
        // cumulative_coverage[6] = fraction at >= 6 = 20 / 100 = 0.2
        assert!((result.cumulative_coverage[6] - 0.2).abs() < 1e-10);
        // cumulative_coverage[10] = fraction at >= 10 = 20 / 100 = 0.2
        assert!((result.cumulative_coverage[10] - 0.2).abs() < 1e-10);
        // cumulative_coverage[11] = fraction at >= 11 = 0 / 100 = 0.0
        assert!((result.cumulative_coverage[11] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_zero_length_chrom() {
        // A chromosome with no bases
        let depths: Vec<u32> = vec![];
        let names = vec!["chrM".to_string()];
        let result = summarize_coverage(&[depths], &names, 100);

        assert_eq!(result.per_chrom.len(), 1);
        let chrom = &result.per_chrom[0];
        assert_eq!(chrom.length, 0);
        assert_eq!(chrom.mean, 0.0);
        assert_eq!(chrom.median, 0.0);
        assert_eq!(chrom.stdev, 0.0);
    }

    #[test]
    fn test_mixed_zero_and_nonzero_chroms() {
        // chr1 is empty, chr2 has 100 bases at depth 5
        let depths1: Vec<u32> = vec![];
        let depths2 = vec![5u32; 100];
        let names = vec!["chr1".to_string(), "chr2".to_string()];
        let result = summarize_coverage(&[depths1, depths2], &names, 100);

        assert_eq!(result.per_chrom.len(), 2);
        assert_eq!(result.per_chrom[0].length, 0);
        assert_eq!(result.per_chrom[1].length, 100);

        // Genome-wide should only reflect chr2
        assert!((result.genome_mean - 5.0).abs() < 1e-10);
        assert!((result.genome_median - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_genome_stdev() {
        // 50 bases at depth 0, 50 bases at depth 10
        // mean = 5.0
        // variance = ((0-5)^2 * 50 + (10-5)^2 * 50) / 100 = (1250 + 1250) / 100 = 25
        // stdev = 5.0
        let mut depths = vec![0u32; 50];
        depths.extend(vec![10u32; 50]);
        let names = vec!["chr1".to_string()];
        let result = summarize_coverage(&[depths], &names, 100);

        assert!((result.genome_stdev - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_max_depth_of_one() {
        // Edge case: max_depth = 1, meaning only one bin (depth 0+)
        let depths = vec![5u32; 10];
        let names = vec!["chr1".to_string()];
        let result = summarize_coverage(&[depths], &names, 1);

        assert_eq!(result.coverage_histogram.len(), 1);
        // All depths capped to bin 0
        assert_eq!(result.coverage_histogram[0], 10);
        assert_eq!(result.cumulative_coverage.len(), 1);
        assert!((result.cumulative_coverage[0] - 1.0).abs() < 1e-10);
    }
}
