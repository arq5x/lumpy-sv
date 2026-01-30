---
name: bedtools-genomecov
description: Compute genome-wide coverage from BED or BAM files using bedtools genomecov.
---

# bedtools genomecov

## Purpose

Compute coverage of a feature file across an entire genome. Produces histograms or per-position depth.

## Key Arguments

| Flag | Description |
|------|-------------|
| `-i` | Input BED/GFF/VCF file |
| `-g` | Genome file (required for BED input). Tab-delimited: chrom, size |
| `-ibam` | Input BAM file (replaces `-i` and `-g`; BAM must be sorted) |
| `-d` | Report depth at each position (1-based). Large output |
| `-dz` | Report depth at each position (0-based), non-zero positions only |
| `-bg` | Output in BedGraph format |
| `-bga` | BedGraph format including zero-coverage regions |
| `-split` | Treat split BAM/BED12 entries as separate intervals |
| `-strand` | Only calculate coverage for `+` or `-` strand |
| `-pc` | Calculate coverage of paired-end fragments (BAM only) |
| `-fs` | Force fragment size (BAM only) |
| `-max` | Cap coverage at this value in histogram |

## Assumptions

- For BED input: genome file (`-g`) is required
- For BAM input: BAM must be position-sorted
- Genome file is tab-delimited: `<chromName>\t<chromSize>`
- Default output is a histogram: chrom, depth, bases_at_depth, chrom_size, fraction

## Common Issues

- **Missing genome file**: BED input requires `-g`. Use `samtools faidx ref.fa` to create a `.fai` file usable as genome file
- **Unsorted BAM**: `-ibam` requires position-sorted BAM. Sort with `samtools sort`
- **Huge output with `-d`**: Reports every position in the genome. Use `-bg` or `-bga` for a more compact representation
- **"genome" line in histogram**: Lines with `genome` in column 1 are the whole-genome summary. Use `grep -w genome` or `grep -vw genome` to isolate

## Validation

Before running:
1. If `-i` is used, confirm genome file exists via `-g`
2. If `-ibam` is used, confirm BAM is sorted (check header with `samtools view -H`)
3. Warn about output size if `-d` is used on a full genome

## Common Patterns

**Genome-wide coverage histogram:**
```bash
bedtools genomecov -i reads.bed -g genome.txt
```

**BedGraph output (good for visualization):**
```bash
bedtools genomecov -ibam reads.bam -bg
```

**BedGraph including zero-coverage regions:**
```bash
bedtools genomecov -ibam reads.bam -bga
```

**Strand-specific coverage:**
```bash
bedtools genomecov -i reads.bed -g genome.txt -strand +
```

## Chaining Hints

- BedGraph output (`-bg`) can be loaded directly into genome browsers (IGV, UCSC)
- `-bga` output piped through `grep -w 0$` extracts all zero-coverage regions
- Often used after alignment to assess sequencing depth
- Output can feed into `bedtools intersect` to find regions above/below coverage thresholds
