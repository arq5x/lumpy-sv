---
name: bedtools-subtract
description: Remove portions of intervals that overlap another file using bedtools subtract.
---

# bedtools subtract

## Purpose

Remove the portion(s) of an interval in A that overlap with features in B. Unlike `intersect -v` which removes entire features, subtract trims them.

## Key Arguments

| Flag | Description |
|------|-------------|
| `-a` | Input file A (features to subtract from) |
| `-b` | Input file B (features to subtract) |
| `-A` | Remove entire A feature if any overlap (instead of trimming) |
| `-N` | Like `-A` but with `-f`, requires sum of all overlaps to meet threshold |
| `-f` | Minimum overlap as fraction of A |
| `-F` | Minimum overlap as fraction of B |
| `-r` | Require reciprocal overlap |
| `-e` | Require overlap fraction for A OR B (not AND) |
| `-s` | Require same strandedness |
| `-S` | Require opposite strandedness |
| `-sorted` | Use sweep algorithm (requires sorted input) |
| `-header` | Print header before results |

## Assumptions

- Both files are tab-delimited with at least 3 columns
- Chromosome names are consistent
- BED coordinates are 0-based, half-open
- If using `-sorted`, both files must be sorted

## Common Issues

- **Subtract vs intersect -v**: `subtract` trims intervals; `intersect -v` removes entire features that have any overlap. Choose based on whether you want trimmed or removed
- **Fragmented output**: Subtracting a feature in the middle of an A feature produces two fragments
- **Using `-A` unexpectedly**: `-A` changes behavior from trimming to complete removal
- **Empty output from chromosome mismatch**: Same issue as intersect

## Validation

Before running:
1. Confirm both files exist
2. Check chromosome naming consistency
3. If `-sorted`, verify both are sorted

## Common Patterns

**Trim away blacklisted regions:**
```bash
bedtools subtract -a peaks.bed -b blacklist.bed
```

**Remove entire features overlapping blacklist:**
```bash
bedtools subtract -a peaks.bed -b blacklist.bed -A
```

**Strand-specific subtraction:**
```bash
bedtools subtract -a genes.bed -b repeats.bed -s
```

## Chaining Hints

- Output may need `bedtools merge` if subtraction creates adjacent fragments
- Often used with blacklist or repeat-masked regions as B
- Alternative to `bedtools intersect -v` when partial removal is desired
