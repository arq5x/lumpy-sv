---
name: bedtools-merge
description: Combine overlapping or nearby intervals into single intervals using bedtools merge.
---

# bedtools merge

## Purpose

Combine overlapping and/or book-ended intervals into a single interval.

## Key Arguments

| Flag | Description |
|------|-------------|
| `-i` | Input file (BED/GFF/VCF). Use `stdin` if piped |
| `-s` | Only merge features on the same strand |
| `-S` | Only merge features on a specific strand (`+` or `-`) |
| `-d` | Max distance between features to merge (default: 0, meaning overlapping and book-ended). Negative values require overlap |
| `-c` | Columns to operate on (comma-delimited list) |
| `-o` | Operations for `-c` columns: sum, min, max, mean, median, mode, collapse, distinct, count, count_distinct, first |
| `-header` | Print header from input before results |
| `-delim` | Delimiter for `collapse` operation (default: comma) |

## Assumptions

- **Input MUST be sorted** by chromosome then start position. This is the most common source of errors
- Tab-delimited with at least 3 columns
- BED coordinates are 0-based, half-open
- If using `-s`, strand column (column 6 in BED) must exist

## Common Issues

- **"ERROR: input is not sorted"**: Merge requires sorted input. Pipe through `bedtools sort` first: `bedtools sort -i file.bed | bedtools merge -i -`
- **Losing metadata**: By default, merge only outputs chrom/start/end. Use `-c` and `-o` to preserve or summarize other columns
- **Book-ended intervals merged**: `chr1 10 20` and `chr1 20 30` will merge into `chr1 10 30` by default. Use `-d -1` to require at least 1bp overlap
- **Unexpected merges**: If `-d` is set too high, distant features will merge. Start with `-d 0` and increase if needed

## Validation

Before running:
1. Confirm input file exists and is sorted. Run: `bedtools sort -i <file> | diff - <file> | head -5`
2. If unsorted, prepend `bedtools sort` to the pipeline
3. If `-c` is used, verify the referenced columns exist

## Common Patterns

**Basic merge of overlapping intervals:**
```bash
bedtools sort -i peaks.bed | bedtools merge -i -
```

**Merge within 100bp and count merged features:**
```bash
bedtools sort -i peaks.bed | bedtools merge -i - -d 100 -c 1 -o count
```

**Merge preserving names (collapsed):**
```bash
bedtools sort -i peaks.bed | bedtools merge -i - -c 4 -o collapse
```

**Strand-specific merge:**
```bash
bedtools sort -i peaks.bed | bedtools merge -i - -s -c 6 -o distinct
```

## Chaining Hints

- Almost always preceded by `bedtools sort`
- Output commonly fed to `bedtools intersect` for downstream overlap analysis
- Use after `bedtools slop` to merge expanded intervals
- Pipe to `bedtools summary` for QC
