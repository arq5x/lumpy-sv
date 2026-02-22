//! Metrics accumulator module.
//!
//! Provides `SampleMetrics`, a struct that accumulates all QC metrics during
//! the single-pass BAM/CRAM iteration over alignment records.

use std::collections::HashMap;

use crate::bam_reader::{CigarOp, ReadRecord};

/// Accumulates all QC metrics for a single sample during a single-pass
/// iteration over a BAM/CRAM file.
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
    pub lengths_for_n50: Vec<u32>,

    // Q-score distribution (mean phred per read -> count)
    pub qscore_histogram: HashMap<u8, u64>,

    // Read length vs quality (subsampled pairs for scatter plot)
    pub length_quality_pairs: Vec<(u32, f32)>,
    pub length_quality_subsample_count: usize,

    // HiFi passes (np tag value -> count)
    pub passes_histogram: HashMap<u32, u64>,
    // Read quality (rq tag, binned: key = (rq * 10000) as u32 -> count)
    pub rq_histogram: HashMap<u32, u64>,

    // Coverage: depth arrays, one per reference sequence
    pub depth_arrays: Vec<Vec<u32>>,

    // Mapping quality distribution (MAPQ -> count)
    pub mapq_histogram: HashMap<u8, u64>,

    // GC content of reads (GC fraction binned to 1%: key = 0..100 -> count)
    pub gc_content_histogram: HashMap<u8, u64>,

    // Mismatch analysis
    pub total_mismatches: u64,

    // Indel analysis
    pub insertion_size_histogram: HashMap<u32, u64>,
    pub deletion_size_histogram: HashMap<u32, u64>,
    pub total_insertions: u64,
    pub total_insertion_bases: u64,
    pub total_deletions: u64,
    pub total_deletion_bases: u64,

    // Sampling flag
    pub is_sampled: bool,
}

impl SampleMetrics {
    /// Create a new `SampleMetrics` accumulator.
    ///
    /// If `skip_coverage` is true, `depth_arrays` will be left empty and no
    /// per-base depth tracking will occur. Otherwise, one `vec![0u32; len]` is
    /// allocated for every reference sequence.
    pub fn new(
        sample_name: String,
        reference_lengths: &[usize],
        skip_coverage: bool,
    ) -> Self {
        let depth_arrays = if skip_coverage {
            Vec::new()
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
            total_mismatches: 0,
            insertion_size_histogram: HashMap::new(),
            deletion_size_histogram: HashMap::new(),
            total_insertions: 0,
            total_insertion_bases: 0,
            total_deletions: 0,
            total_deletion_bases: 0,
            is_sampled: false,
        }
    }

    /// Process a single read, updating every accumulator.
    pub fn process_read(&mut self, read: &ReadRecord) {
        // Basic counts
        self.total_reads += 1;
        let seq_len = read.sequence_length as u32;
        self.total_bases += seq_len as u64;

        if read.is_unmapped {
            self.unmapped_reads += 1;
        } else {
            self.mapped_reads += 1;
        }

        if read.is_duplicate {
            self.duplicate_reads += 1;
        }

        // Read length distribution
        *self.length_histogram.entry(seq_len).or_insert(0) += 1;
        self.lengths_for_n50.push(seq_len);

        // Mean Q-score for this read
        let mean_q = Self::compute_mean_quality(&read.quality_scores);
        *self.qscore_histogram.entry(mean_q).or_insert(0) += 1;

        // Subsample length-quality pairs (keep first N)
        if self.length_quality_pairs.len() < self.length_quality_subsample_count {
            self.length_quality_pairs.push((seq_len, mean_q as f32));
        }

        // Mapping quality
        *self.mapq_histogram.entry(read.mapping_quality).or_insert(0) += 1;

        // HiFi passes (np tag)
        if let Some(np) = read.num_passes {
            if np >= 0 {
                *self.passes_histogram.entry(np as u32).or_insert(0) += 1;
            }
        }

        // Read quality (rq tag, binned)
        if let Some(rq) = read.read_quality {
            let binned = (rq * 10_000.0) as u32;
            *self.rq_histogram.entry(binned).or_insert(0) += 1;
        }

        // GC content
        if !read.sequence_bases.is_empty() {
            let gc_count = read
                .sequence_bases
                .iter()
                .filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
                .count();
            let gc_fraction = gc_count as f64 / read.sequence_bases.len() as f64;
            let gc_bin = (gc_fraction * 100.0).round().min(100.0) as u8;
            *self.gc_content_histogram.entry(gc_bin).or_insert(0) += 1;
        }

        // Alignment-derived metrics (mapped reads only)
        if !read.is_unmapped {
            if let (Some(ref_id), Some(start)) = (read.reference_id, read.alignment_start) {
                self.process_alignment(ref_id, start, &read.cigar_ops);
            }
        }
    }

    /// Walk the CIGAR string to update depth arrays, mismatch counts, and
    /// indel statistics.
    fn process_alignment(
        &mut self,
        ref_id: usize,
        start: usize,
        cigar_ops: &[(CigarOp, usize)],
    ) {
        let mut ref_pos = start;

        for &(op, len) in cigar_ops {
            match op {
                // Reference-consuming + depth-contributing ops
                CigarOp::Match | CigarOp::SequenceMatch => {
                    self.increment_depth(ref_id, ref_pos, len);
                    ref_pos += len;
                }
                CigarOp::SequenceMismatch => {
                    self.increment_depth(ref_id, ref_pos, len);
                    self.total_mismatches += len as u64;
                    ref_pos += len;
                }
                // Deletions consume reference but do not contribute depth
                CigarOp::Deletion => {
                    self.total_deletions += 1;
                    self.total_deletion_bases += len as u64;
                    *self
                        .deletion_size_histogram
                        .entry(len as u32)
                        .or_insert(0) += 1;
                    ref_pos += len;
                }
                // Insertions consume query but not reference
                CigarOp::Insertion => {
                    self.total_insertions += 1;
                    self.total_insertion_bases += len as u64;
                    *self
                        .insertion_size_histogram
                        .entry(len as u32)
                        .or_insert(0) += 1;
                }
                // Skip (N) consumes reference
                CigarOp::Skip => {
                    ref_pos += len;
                }
                // Soft clip, hard clip, and pad do not consume reference
                CigarOp::SoftClip | CigarOp::HardClip | CigarOp::Pad => {}
            }
        }
    }

    /// Increment depth values over `[start, start + len)` for a given
    /// reference sequence. Uses `saturating_add` to avoid overflow.
    fn increment_depth(&mut self, ref_id: usize, start: usize, len: usize) {
        if let Some(depth) = self.depth_arrays.get_mut(ref_id) {
            let end = (start + len).min(depth.len());
            let begin = start.min(depth.len());
            for pos in begin..end {
                depth[pos] = depth[pos].saturating_add(1);
            }
        }
    }

    /// Compute the mean Phred quality score for a set of quality values,
    /// returned as a rounded `u8`.
    fn compute_mean_quality(quality_scores: &[u8]) -> u8 {
        if quality_scores.is_empty() {
            return 0;
        }
        let sum: u64 = quality_scores.iter().map(|&q| q as u64).sum();
        let mean = sum as f64 / quality_scores.len() as f64;
        mean.round() as u8
    }

    /// Compute N50: the read length such that reads of this length or longer
    /// account for at least 50% of total bases.
    pub fn compute_n50(&self) -> u32 {
        if self.lengths_for_n50.is_empty() {
            return 0;
        }

        let mut sorted = self.lengths_for_n50.clone();
        sorted.sort_unstable_by(|a, b| b.cmp(a)); // descending

        let total: u64 = sorted.iter().map(|&l| l as u64).sum();
        let half = total / 2;

        let mut cumulative: u64 = 0;
        for &length in &sorted {
            cumulative += length as u64;
            if cumulative >= half {
                return length;
            }
        }

        // Should not reach here, but return the smallest length as fallback
        *sorted.last().unwrap_or(&0)
    }

    /// Mean read length across all processed reads.
    pub fn mean_read_length(&self) -> f64 {
        if self.total_reads == 0 {
            return 0.0;
        }
        self.total_bases as f64 / self.total_reads as f64
    }

    /// Median read length across all processed reads.
    pub fn median_read_length(&self) -> u32 {
        if self.lengths_for_n50.is_empty() {
            return 0;
        }

        let mut sorted = self.lengths_for_n50.clone();
        sorted.sort_unstable();

        let mid = sorted.len() / 2;
        if sorted.len() % 2 == 0 {
            // Average of the two middle values
            ((sorted[mid - 1] as u64 + sorted[mid] as u64) / 2) as u32
        } else {
            sorted[mid]
        }
    }

    /// Overall mean Q-score across all reads, weighted by read count per
    /// quality bin.
    pub fn mean_qscore(&self) -> f64 {
        let total_count: u64 = self.qscore_histogram.values().sum();
        if total_count == 0 {
            return 0.0;
        }
        let weighted_sum: u64 = self
            .qscore_histogram
            .iter()
            .map(|(&q, &count)| q as u64 * count)
            .sum();
        weighted_sum as f64 / total_count as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam_reader::{CigarOp, ReadRecord};

    /// Helper to build a minimal ReadRecord for testing.
    fn make_read(
        seq_len: usize,
        quals: Vec<u8>,
        mapped: bool,
        cigar: Vec<(CigarOp, usize)>,
    ) -> ReadRecord {
        let bases = vec![b'A'; seq_len];
        ReadRecord {
            sequence_length: seq_len,
            quality_scores: quals,
            mapping_quality: 60,
            reference_id: if mapped { Some(0) } else { None },
            alignment_start: if mapped { Some(0) } else { None },
            cigar_ops: cigar,
            sequence_bases: bases,
            num_passes: Some(10),
            read_quality: Some(0.999),
            is_unmapped: !mapped,
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        }
    }

    #[test]
    fn test_new_with_coverage() {
        let m = SampleMetrics::new("sample1".to_string(), &[100, 200], false);
        assert_eq!(m.depth_arrays.len(), 2);
        assert_eq!(m.depth_arrays[0].len(), 100);
        assert_eq!(m.depth_arrays[1].len(), 200);
        assert_eq!(m.sample_name, "sample1");
    }

    #[test]
    fn test_new_skip_coverage() {
        let m = SampleMetrics::new("s".to_string(), &[100, 200], true);
        assert!(m.depth_arrays.is_empty());
    }

    #[test]
    fn test_process_unmapped_read() {
        let mut m = SampleMetrics::new("s".to_string(), &[], true);
        let read = make_read(1000, vec![30; 1000], false, vec![]);
        m.process_read(&read);

        assert_eq!(m.total_reads, 1);
        assert_eq!(m.unmapped_reads, 1);
        assert_eq!(m.mapped_reads, 0);
        assert_eq!(m.total_bases, 1000);
        assert_eq!(*m.length_histogram.get(&1000).unwrap(), 1);
    }

    #[test]
    fn test_process_mapped_read_with_depth() {
        let mut m = SampleMetrics::new("s".to_string(), &[1000], false);
        let cigar = vec![(CigarOp::SequenceMatch, 500)];
        let read = make_read(500, vec![20; 500], true, cigar);
        m.process_read(&read);

        assert_eq!(m.mapped_reads, 1);
        // Check depth at position 0 and 499
        assert_eq!(m.depth_arrays[0][0], 1);
        assert_eq!(m.depth_arrays[0][499], 1);
        assert_eq!(m.depth_arrays[0][500], 0);
    }

    #[test]
    fn test_indel_tracking() {
        let mut m = SampleMetrics::new("s".to_string(), &[1000], false);
        let cigar = vec![
            (CigarOp::SequenceMatch, 100),
            (CigarOp::Insertion, 5),
            (CigarOp::SequenceMatch, 50),
            (CigarOp::Deletion, 3),
            (CigarOp::SequenceMatch, 50),
        ];
        let read = make_read(205, vec![30; 205], true, cigar);
        m.process_read(&read);

        assert_eq!(m.total_insertions, 1);
        assert_eq!(m.total_insertion_bases, 5);
        assert_eq!(m.total_deletions, 1);
        assert_eq!(m.total_deletion_bases, 3);
        assert_eq!(*m.insertion_size_histogram.get(&5).unwrap(), 1);
        assert_eq!(*m.deletion_size_histogram.get(&3).unwrap(), 1);
    }

    #[test]
    fn test_mismatch_tracking() {
        let mut m = SampleMetrics::new("s".to_string(), &[1000], false);
        let cigar = vec![
            (CigarOp::SequenceMatch, 90),
            (CigarOp::SequenceMismatch, 10),
        ];
        let read = make_read(100, vec![30; 100], true, cigar);
        m.process_read(&read);

        assert_eq!(m.total_mismatches, 10);
    }

    #[test]
    fn test_compute_n50() {
        let mut m = SampleMetrics::new("s".to_string(), &[], true);
        // 5 reads: lengths 100, 200, 300, 400, 500
        // total = 1500, half = 750
        // sorted desc: 500, 400, 300, 200, 100
        // cumulative: 500 < 750, 500+400=900 >= 750 -> N50 = 400
        for len in [100u32, 200, 300, 400, 500] {
            m.lengths_for_n50.push(len);
        }
        assert_eq!(m.compute_n50(), 400);
    }

    #[test]
    fn test_mean_read_length() {
        let mut m = SampleMetrics::new("s".to_string(), &[], true);
        m.total_reads = 4;
        m.total_bases = 1000;
        assert!((m.mean_read_length() - 250.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_median_read_length_odd() {
        let mut m = SampleMetrics::new("s".to_string(), &[], true);
        m.lengths_for_n50 = vec![100, 300, 200];
        assert_eq!(m.median_read_length(), 200);
    }

    #[test]
    fn test_median_read_length_even() {
        let mut m = SampleMetrics::new("s".to_string(), &[], true);
        m.lengths_for_n50 = vec![100, 200, 300, 400];
        // median = (200 + 300) / 2 = 250
        assert_eq!(m.median_read_length(), 250);
    }

    #[test]
    fn test_mean_qscore() {
        let mut m = SampleMetrics::new("s".to_string(), &[], true);
        m.qscore_histogram.insert(20, 5);
        m.qscore_histogram.insert(30, 5);
        // (20*5 + 30*5) / 10 = 250/10 = 25.0
        assert!((m.mean_qscore() - 25.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_gc_content() {
        let mut m = SampleMetrics::new("s".to_string(), &[], true);
        // Build a read with 50% GC content
        let mut bases = vec![b'G'; 50];
        bases.extend(vec![b'A'; 50]);
        let read = ReadRecord {
            sequence_length: 100,
            quality_scores: vec![30; 100],
            mapping_quality: 60,
            reference_id: None,
            alignment_start: None,
            cigar_ops: vec![],
            sequence_bases: bases,
            num_passes: None,
            read_quality: None,
            is_unmapped: true,
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        m.process_read(&read);
        assert_eq!(*m.gc_content_histogram.get(&50).unwrap(), 1);
    }

    #[test]
    fn test_length_quality_subsample() {
        let mut m = SampleMetrics::new("s".to_string(), &[], true);
        m.length_quality_subsample_count = 5;
        for _ in 0..10 {
            let read = make_read(100, vec![30; 100], false, vec![]);
            m.process_read(&read);
        }
        assert_eq!(m.length_quality_pairs.len(), 5);
    }

    #[test]
    fn test_empty_metrics() {
        let m = SampleMetrics::new("empty".to_string(), &[], true);
        assert_eq!(m.compute_n50(), 0);
        assert_eq!(m.mean_read_length(), 0.0);
        assert_eq!(m.median_read_length(), 0);
        assert_eq!(m.mean_qscore(), 0.0);
    }

    #[test]
    fn test_depth_saturation() {
        let mut m = SampleMetrics::new("s".to_string(), &[10], false);
        // Set depth to near u32::MAX
        m.depth_arrays[0][0] = u32::MAX - 1;
        // Increment depth at position 0 twice: should saturate at u32::MAX
        m.increment_depth(0, 0, 1);
        assert_eq!(m.depth_arrays[0][0], u32::MAX);
        m.increment_depth(0, 0, 1);
        assert_eq!(m.depth_arrays[0][0], u32::MAX); // saturated
    }

    #[test]
    fn test_duplicate_tracking() {
        let mut m = SampleMetrics::new("s".to_string(), &[], true);
        let mut read = make_read(100, vec![30; 100], true, vec![]);
        read.is_duplicate = true;
        m.process_read(&read);
        assert_eq!(m.duplicate_reads, 1);
    }

    #[test]
    fn test_passes_and_rq_histograms() {
        let mut m = SampleMetrics::new("s".to_string(), &[], true);
        let read = make_read(100, vec![30; 100], false, vec![]);
        m.process_read(&read);

        // np tag was 10
        assert_eq!(*m.passes_histogram.get(&10).unwrap(), 1);
        // rq tag was 0.999 -> 0.999 * 10000 = 9990
        assert_eq!(*m.rq_histogram.get(&9990).unwrap(), 1);
    }
}
