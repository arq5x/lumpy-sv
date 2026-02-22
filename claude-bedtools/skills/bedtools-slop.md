---
name: bedtools-slop
description: Extend intervals by a specified number of bases using bedtools slop.
---

# bedtools slop

## Purpose

Add "slop" (padding) to each feature — extend intervals upstream, downstream, or both directions.

## Key Arguments

| Flag | Description |
|------|-------------|
| `-i` | Input file (BED/GFF/VCF) |
| `-g` | Genome file (required). Prevents coordinates from exceeding chromosome bounds |
| `-b` | Extend both directions by this many bases |
| `-l` | Extend left (subtract from start) by this many bases |
| `-r` | Extend right (add to end) by this many bases |
| `-s` | Use strand to define left/right (left=upstream, right=downstream) |
| `-pct` | Treat `-l`/`-r`/`-b` as fraction of feature length |
| `-header` | Print header before results |

## Assumptions

- Genome file is required to prevent coordinates from going out of bounds
- Tab-delimited with at least 3 columns
- Starts are clipped to 0; ends are clipped to chromosome length
- If using `-s`, strand column (column 6 in BED) must exist
- Either `-b` or both `-l` and `-r` must be specified

## Common Issues

- **Missing genome file**: Always required. Create from FASTA index: `samtools faidx ref.fa`
- **Confusing `-l` and `-r` with strand**: Without `-s`, `-l` always subtracts from start and `-r` adds to end. With `-s`, they become upstream/downstream relative to strand
- **Overlapping results**: Slop can cause previously non-overlapping intervals to overlap. Pipe to `bedtools merge` if needed
- **Using `-pct` with `-b`**: `-pct` interprets the value as a fraction. `-b 0.1 -pct` on a 1000bp feature adds 100bp each side

## Validation

Before running:
1. Confirm input and genome files exist
2. Verify either `-b` or both `-l` and `-r` are specified
3. If `-s` is used, verify strand column exists

## Common Patterns

**Add 500bp flanks both sides:**
```bash
bedtools slop -i peaks.bed -g genome.txt -b 500
```

**Add 2kb upstream of genes (strand-aware):**
```bash
bedtools slop -i genes.bed -g genome.txt -l 2000 -r 0 -s
```

**Extend by 10% of feature length:**
```bash
bedtools slop -i features.bed -g genome.txt -b 0.1 -pct
```

## Chaining Hints

- Output often piped to `bedtools merge` to combine newly overlapping intervals
- Used before `bedtools getfasta` to extract flanking sequences
- Used before `bedtools intersect` to find features within a distance of targets
- Alternative to `bedtools closest` for "within X bp" queries: slop + intersect
