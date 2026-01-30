---
name: bedtools-sort
description: Sort genomic intervals by chromosome and position using bedtools sort.
---

# bedtools sort

## Purpose

Sort a feature file by chromosome and start position (default), or by other criteria like feature size or score.

## Key Arguments

| Flag | Description |
|------|-------------|
| `-i` | Input file (BED/GFF/VCF) |
| `-sizeA` | Sort by feature size, ascending |
| `-sizeD` | Sort by feature size, descending |
| `-chrThenSizeA` | Sort by chrom (asc), then size (asc) |
| `-chrThenSizeD` | Sort by chrom (asc), then size (desc) |
| `-chrThenScoreA` | Sort by chrom (asc), then score (asc) |
| `-chrThenScoreD` | Sort by chrom (asc), then score (desc) |
| `-g` | Sort by chromosome order in a genome file |
| `-faidx` | Sort by chromosome order in a FASTA index file |
| `-header` | Print header line(s) before results |

## Assumptions

- Tab-delimited with at least 3 columns
- If using `-g` or `-faidx`, the genome/index file must exist and list all chromosomes present in the input
- If using score-based sorting, column 5 must contain numeric scores

## Common Issues

- **Default sort is lexicographic for chromosomes**: `chr10` sorts before `chr2`. Use `-g` or `-faidx` for natural chromosome order
- **Header lines get sorted into the body**: Use `-header` to preserve header lines at the top
- **Sort order differs from `sort -k1,1 -k2,2n`**: bedtools sort and Unix sort may produce different chromosome orders. For bedtools operations, always use `bedtools sort`

## Validation

Before running:
1. Confirm input file exists
2. If `-g` or `-faidx` is specified, confirm the file exists

## Common Patterns

**Default sort (chrom + start):**
```bash
bedtools sort -i unsorted.bed
```

**Sort by genome file chromosome order:**
```bash
bedtools sort -i unsorted.bed -g genome.txt
```

**Sort by feature size:**
```bash
bedtools sort -i features.bed -sizeD
```

## Chaining Hints

- Required before `bedtools merge`, `bedtools complement`, and `bedtools closest`
- Required before using `-sorted` flag in `bedtools intersect`
- Commonly the first step in any pipeline
