---
name: bedtools-complement
description: Find regions of the genome NOT covered by intervals using bedtools complement.
---

# bedtools complement

## Purpose

Return the base pair complement of a feature file — i.e., all genomic regions NOT covered by the input intervals.

## Key Arguments

| Flag | Description |
|------|-------------|
| `-i` | Input file (BED/GFF/VCF) |
| `-g` | Genome file (required). Tab-delimited: chrom, size |
| `-L` | Limit output to chromosomes present in the input file |

## Assumptions

- **Input MUST be sorted** by chromosome then start position
- Genome file is required and must list all chromosomes and their sizes
- Tab-delimited input with at least 3 columns
- Genome file format: `<chromName>\t<chromSize>`

## Common Issues

- **Missing genome file**: Always required. Create from FASTA index: `samtools faidx ref.fa` then use the `.fai` file
- **Unsorted input**: Will produce errors. Always sort first
- **Chromosomes in input not in genome file**: Will cause errors. Ensure genome file covers all chromosomes in input
- **Overlapping input intervals**: Complement works on the merged footprint, but unsorted overlapping intervals may cause errors. Sort and merge first if input has overlaps

## Validation

Before running:
1. Confirm input file exists and is sorted
2. Confirm genome file exists
3. Check that all chromosomes in input appear in genome file

## Common Patterns

**Find uncovered regions:**
```bash
bedtools sort -i features.bed | bedtools complement -i - -g genome.txt
```

**Complement limited to chromosomes with data:**
```bash
bedtools sort -i features.bed | bedtools complement -i - -g genome.txt -L
```

## Chaining Hints

- Always preceded by `bedtools sort` (and optionally `bedtools merge`)
- Output is useful as input to `bedtools getfasta` to extract sequences of uncovered regions
- Complement of exons = introns + intergenic regions
- Complement of peaks = non-peak regions for background analysis
