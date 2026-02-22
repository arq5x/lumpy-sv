---
name: bedtools-summary
description: Report summary statistics and QC for interval files using bedtools summary.
---

# bedtools summary

## Purpose

Report summary statistics of the intervals in a file. Useful for quality control and quick characterization of genomic interval data.

## Key Arguments

| Flag | Description |
|------|-------------|
| `-i` | Input file (BED/GFF/VCF/BAM) |
| `-g` | Genome file (required). Tab-delimited: chrom, size |

## Assumptions

- Genome file is required
- Input is tab-delimited with at least 3 columns (for BED/GFF/VCF)
- Genome file format: `<chromName>\t<chromSize>`

## Common Issues

- **Missing genome file**: Always required. Create from FASTA index: `samtools faidx ref.fa`
- **Chromosomes in input not in genome file**: May cause errors or incomplete summary

## Validation

Before running:
1. Confirm input file exists
2. Confirm genome file exists
3. Check that chromosomes in input appear in genome file

## Common Patterns

**Basic summary of a BED file:**
```bash
bedtools summary -i features.bed -g genome.txt
```

**QC after a pipeline:**
```bash
bedtools intersect -a peaks.bed -b genes.bed > result.bed
bedtools summary -i result.bed -g genome.txt
```

## Chaining Hints

- Use as the final step in any pipeline for QC
- Suggested after complex multi-step operations to verify results are reasonable
- Run on both input and output files to compare before/after statistics
