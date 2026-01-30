---
name: bedtools-closest
description: Find the nearest feature in another file for each interval using bedtools closest.
---

# bedtools closest

## Purpose

For each feature in A, find the closest feature (upstream or downstream) in B.

## Key Arguments

| Flag | Description |
|------|-------------|
| `-a` | Input file A (query) |
| `-b` | Input file B (database) |
| `-d` | Report distance to closest feature as extra column (0 for overlaps) |
| `-D` | Like `-d` but with signed distance. Options: `ref`, `a`, `b` for orientation |
| `-io` | Ignore overlapping features (only report nearby, non-touching) |
| `-iu` | Ignore upstream features (requires `-D`) |
| `-id` | Ignore downstream features (requires `-D`) |
| `-fu` | Choose first upstream feature (requires `-D`) |
| `-fd` | Choose first downstream feature (requires `-D`) |
| `-t` | Tie handling: `all` (default), `first`, `last` |
| `-s` | Require same strandedness |
| `-S` | Require opposite strandedness |
| `-k` | Report the k closest features (default: 1) |
| `-header` | Print header from A before results |

## Assumptions

- **Both files MUST be sorted** by chromosome then start position
- Tab-delimited with at least 3 columns
- Chromosome names must be consistent across files
- If using `-D` with `a` or `b` orientation, strand column must exist

## Common Issues

- **Unsorted input**: closest requires sorted input and will produce wrong results silently if not sorted
- **Features on different chromosomes**: If a feature in A has no match on the same chromosome in B, the output will show `.` and distance `-1`
- **Tie handling**: When multiple B features are equidistant, all are reported by default. Use `-t first` to get one result per A feature
- **Signed vs unsigned distance**: `-d` gives absolute distance; `-D` gives signed distance where negative means upstream

## Validation

Before running:
1. Confirm both files exist and are sorted
2. Check chromosome naming consistency
3. If `-D` is used with orientation, verify strand column exists

## Common Patterns

**Find nearest gene for each peak:**
```bash
bedtools sort -i peaks.bed > sorted_peaks.bed
bedtools sort -i genes.bed > sorted_genes.bed
bedtools closest -a sorted_peaks.bed -b sorted_genes.bed -d
```

**Find nearest non-overlapping upstream feature:**
```bash
bedtools closest -a sorted_a.bed -b sorted_b.bed -D ref -iu -io
```

**Find 3 closest features:**
```bash
bedtools closest -a sorted_a.bed -b sorted_b.bed -k 3 -d
```

## Chaining Hints

- Always preceded by `bedtools sort` for both inputs
- Often used after `bedtools intersect -v` to find nearest feature for non-overlapping entries
- Pipe to `awk` or `bedtools summary` to analyze distance distributions
