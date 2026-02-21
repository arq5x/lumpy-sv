//! BAM/CRAM reader module.
//!
//! Provides a unified interface for reading BAM and CRAM files using the noodles crate.
//! Extracts header information and iterates over alignment records, converting them
//! into our own `ReadRecord` type for downstream QC processing.

use std::fs::File;
use std::io;
use std::path::Path;

use noodles::bam;
use noodles::cram;
use noodles::fasta;
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// CIGAR operation kinds, mirroring noodles' `Kind` enum.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match,
    Insertion,
    Deletion,
    Skip,
    SoftClip,
    HardClip,
    Pad,
    SequenceMatch,
    SequenceMismatch,
}

impl From<Kind> for CigarOp {
    fn from(kind: Kind) -> Self {
        match kind {
            Kind::Match => CigarOp::Match,
            Kind::Insertion => CigarOp::Insertion,
            Kind::Deletion => CigarOp::Deletion,
            Kind::Skip => CigarOp::Skip,
            Kind::SoftClip => CigarOp::SoftClip,
            Kind::HardClip => CigarOp::HardClip,
            Kind::Pad => CigarOp::Pad,
            Kind::SequenceMatch => CigarOp::SequenceMatch,
            Kind::SequenceMismatch => CigarOp::SequenceMismatch,
        }
    }
}

/// Information extracted from the BAM/CRAM header.
#[derive(Clone, Debug)]
pub struct FileHeader {
    /// Chromosome / reference sequence names.
    pub reference_names: Vec<String>,
    /// Chromosome / reference sequence lengths.
    pub reference_lengths: Vec<usize>,
    /// Sample name from the first read group's SM tag, if present.
    pub sample_name: Option<String>,
}

/// Per-read data extracted during iteration.
#[derive(Clone, Debug)]
pub struct ReadRecord {
    pub sequence_length: usize,
    pub quality_scores: Vec<u8>,
    pub mapping_quality: u8,
    pub reference_id: Option<usize>,
    /// 0-based alignment start position.
    pub alignment_start: Option<usize>,
    pub cigar_ops: Vec<(CigarOp, usize)>,
    /// Sequence bases (A, C, G, T, N, etc.).
    pub sequence_bases: Vec<u8>,
    /// Number of passes (PacBio `np` tag).
    pub num_passes: Option<i64>,
    /// Read quality (PacBio `rq` tag).
    pub read_quality: Option<f32>,
    pub is_unmapped: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
    pub is_duplicate: bool,
}

// ---------------------------------------------------------------------------
// Header extraction
// ---------------------------------------------------------------------------

/// Extract a `FileHeader` from a noodles `sam::Header`.
fn extract_header(header: &sam::Header) -> FileHeader {
    let mut reference_names = Vec::new();
    let mut reference_lengths = Vec::new();

    for (name, ref_seq) in header.reference_sequences() {
        let name_str = std::str::from_utf8(name.as_ref())
            .unwrap_or("unknown")
            .to_string();
        reference_names.push(name_str);
        reference_lengths.push(usize::from(ref_seq.length()));
    }

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

/// Read just the header from a BAM or CRAM file.
///
/// The `reference` parameter is optional; it is used for CRAM files
/// that need a reference FASTA to decode records.
pub fn read_header<P: AsRef<Path>>(
    path: P,
    reference: Option<&str>,
) -> io::Result<FileHeader> {
    let path = path.as_ref();
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("");

    match ext {
        "cram" => {
            let repo = build_fasta_repo(reference)?;
            let mut reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(path)?;
            let header = reader.read_header()?;
            Ok(extract_header(&header))
        }
        _ => {
            // Treat everything else (bam, etc.) as BAM
            let mut reader = File::open(path).map(bam::io::Reader::new)?;
            let header = reader.read_header()?;
            Ok(extract_header(&header))
        }
    }
}

// ---------------------------------------------------------------------------
// Record extraction helpers
// ---------------------------------------------------------------------------

/// Tag constants for PacBio-specific tags.
const TAG_NP: Tag = Tag::new(b'n', b'p');
const TAG_RQ: Tag = Tag::new(b'r', b'q');

/// Extract a `ReadRecord` from a BAM record.
fn read_record_from_bam(record: &bam::Record) -> io::Result<ReadRecord> {
    let flags = record.flags();

    let is_unmapped = flags.is_unmapped();
    let is_secondary = flags.is_secondary();
    let is_supplementary = flags.is_supplementary();
    let is_duplicate = flags.is_duplicate();

    let mapping_quality = record
        .mapping_quality()
        .map(|mq| u8::from(mq))
        .unwrap_or(255);

    let reference_id = record
        .reference_sequence_id()
        .map(|r| r)
        .transpose()?;

    let alignment_start = record
        .alignment_start()
        .map(|r| r)
        .transpose()?
        .map(|pos| usize::from(pos) - 1); // Convert 1-based to 0-based

    // Sequence
    let seq = record.sequence();
    let sequence_length = seq.len();
    let sequence_bases: Vec<u8> = seq.iter().collect();

    // Quality scores (raw Phred, not ASCII-offset)
    let qual = record.quality_scores();
    let quality_scores: Vec<u8> = qual.as_ref().to_vec();

    // CIGAR
    let cigar = record.cigar();
    let mut cigar_ops = Vec::with_capacity(cigar.len());
    for op_result in cigar.iter() {
        let op = op_result?;
        cigar_ops.push((CigarOp::from(op.kind()), op.len()));
    }

    // Tags
    let data = record.data();
    let num_passes = extract_int_tag_bam(&data, &TAG_NP);
    let read_quality = extract_float_tag_bam(&data, &TAG_RQ);

    Ok(ReadRecord {
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
    })
}

/// Extract an integer tag value from BAM record data.
fn extract_int_tag_bam(data: &bam::record::Data<'_>, tag: &Tag) -> Option<i64> {
    data.get(tag)
        .and_then(|result| result.ok())
        .and_then(|value| value.as_int())
}

/// Extract a float tag value from BAM record data.
fn extract_float_tag_bam(
    data: &bam::record::Data<'_>,
    tag: &Tag,
) -> Option<f32> {
    use noodles::sam::alignment::record::data::field::Value;

    data.get(tag)
        .and_then(|result| result.ok())
        .and_then(|value| match value {
            Value::Float(f) => Some(f),
            _ => value.as_int().map(|i| i as f32),
        })
}

/// Extract a `ReadRecord` from a CRAM record.
fn read_record_from_cram(
    record: &cram::Record,
    header: &sam::Header,
) -> io::Result<ReadRecord> {
    let flags = record.flags();

    let is_unmapped = flags.is_unmapped();
    let is_secondary = flags.is_secondary();
    let is_supplementary = flags.is_supplementary();
    let is_duplicate = flags.is_duplicate();

    let mapping_quality = record
        .mapping_quality()
        .map(|mq| u8::from(mq))
        .unwrap_or(255);

    let reference_id = record.reference_sequence_id();

    let alignment_start = record
        .alignment_start()
        .map(|pos| usize::from(pos) - 1); // Convert 1-based to 0-based

    // Sequence
    let seq = record.sequence();
    let sequence_length = seq.len();
    let sequence_bases: Vec<u8> = seq.as_ref().to_vec();

    // Quality scores (raw Phred)
    let qual = record.quality_scores();
    let quality_scores: Vec<u8> = qual.as_ref().to_vec();

    // CIGAR - use the sam::alignment::Record trait implementation
    let cigar_trait: Box<dyn noodles::sam::alignment::record::Cigar + '_> =
        sam::alignment::Record::cigar(record);
    let mut cigar_ops = Vec::new();
    for op_result in cigar_trait.iter() {
        let op = op_result?;
        cigar_ops.push((CigarOp::from(op.kind()), op.len()));
    }

    // Tags (CRAM uses record_buf::Data)
    let data = record.data();
    let num_passes = extract_int_tag_cram(data, &TAG_NP);
    let read_quality = extract_float_tag_cram(data, &TAG_RQ);

    // Suppress unused variable warning
    let _ = header;

    Ok(ReadRecord {
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
    })
}

/// Extract an integer tag value from CRAM record data (record_buf::Data).
fn extract_int_tag_cram(
    data: &sam::alignment::record_buf::Data,
    tag: &Tag,
) -> Option<i64> {
    data.get(tag).and_then(|value| value.as_int())
}

/// Extract a float tag value from CRAM record data (record_buf::Data).
fn extract_float_tag_cram(
    data: &sam::alignment::record_buf::Data,
    tag: &Tag,
) -> Option<f32> {
    use noodles::sam::alignment::record_buf::data::field::Value;

    data.get(tag).and_then(|value| match value {
        Value::Float(f) => Some(*f),
        _ => value.as_int().map(|i| i as f32),
    })
}

// ---------------------------------------------------------------------------
// FASTA reference repository builder
// ---------------------------------------------------------------------------

/// Build a FASTA reference sequence repository from an optional path.
fn build_fasta_repo(reference: Option<&str>) -> io::Result<fasta::Repository> {
    match reference {
        Some(ref_path) => {
            let indexed_reader = fasta::io::indexed_reader::Builder::default()
                .build_from_path(ref_path)?;
            let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
            Ok(fasta::Repository::new(adapter))
        }
        None => Ok(fasta::Repository::default()),
    }
}

// ---------------------------------------------------------------------------
// Main processing entry points
// ---------------------------------------------------------------------------

/// Open a BAM or CRAM file and iterate over all PRIMARY (non-secondary,
/// non-supplementary) records, calling `callback` for each one.
///
/// The file format is determined by the file extension:
/// - `.cram` -> CRAM reader (with optional reference)
/// - anything else -> BAM reader
pub fn process_alignment_file<P, F>(
    path: P,
    reference: Option<&str>,
    mut callback: F,
) -> io::Result<()>
where
    P: AsRef<Path>,
    F: FnMut(ReadRecord),
{
    let path = path.as_ref();
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("");

    match ext {
        "cram" => process_cram(path, reference, &mut callback),
        _ => process_bam(path, &mut callback),
    }
}

/// Process a BAM file, calling `callback` for each primary record.
fn process_bam<F>(path: &Path, callback: &mut F) -> io::Result<()>
where
    F: FnMut(ReadRecord),
{
    let mut reader = File::open(path).map(bam::io::Reader::new)?;
    let _header = reader.read_header()?;

    for result in reader.records() {
        let record = result?;
        let flags = record.flags();

        // Skip secondary and supplementary alignments
        if flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        let read_record = read_record_from_bam(&record)?;
        callback(read_record);
    }

    Ok(())
}

/// Process a CRAM file, calling `callback` for each primary record.
fn process_cram<F>(
    path: &Path,
    reference: Option<&str>,
    callback: &mut F,
) -> io::Result<()>
where
    F: FnMut(ReadRecord),
{
    let repo = build_fasta_repo(reference)?;
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repo)
        .build_from_path(path)?;
    let header = reader.read_header()?;

    for result in reader.records(&header) {
        let record = result?;
        let flags = record.flags();

        // Skip secondary and supplementary alignments
        if flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        let read_record = read_record_from_cram(&record, &header)?;
        callback(read_record);
    }

    Ok(())
}
