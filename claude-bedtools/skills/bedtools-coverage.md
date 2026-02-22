---
name: bedtools-coverage
description: Compute depth and breadth of coverage of B features on A intervals using bedtools coverage.
---

# bedtools coverage

## Purpose

Compute the depth and breadth of coverage of features from B on the intervals in A.

## Key Arguments

| Flag | Description |
|------|-------------|
| `-a` | Input file A (intervals to measure coverage on) |
| `-b` | Input file B (features providing coverage) |
| `-hist` | Report histogram of coverage per A feature and summary |
| `-d` | Report depth at each position in each A feature (1-based) |
| `-counts` | Only report overlap count, skip fraction computation |
| `-mean` | Report mean depth across each A feature |
| `-s` | Require same strandedness |
| `-S` | Require opposite strandedness |
| `-f` | Minimum overlap as fraction of A |
| `-F` | Minimum overlap as fraction of B |
| `-r` | Require reciprocal overlap |
| `-sorted` | Use sweep algorithm (requires sorted input) |
| `-header` | Print header from A before results |

## Assumptions

- Both files are tab-delimited with at least 3 columns
- Chromosome names are consistent
- If using `-sorted`, both files must be sorted
- Default output appends 4 columns to each A entry: count of overlaps, number of bases covered, length of A, fraction of A covered

## Common Issues

- **Confusing A and B**: A is the "windows" you measure coverage on; B provides the "reads" or features. Swapping them gives very different results
- **Output format varies by flag**: Default, `-hist`, `-d`, `-counts`, and `-mean` all produce different output formats
- **Large output with `-d`**: Reports one line per base per A feature. For large features, output can be enormous
- **Memory with unsorted input**: Use `-sorted` for large files

## Validation

Before running:
1. Confirm both files exist
2. If `-sorted`, verify both are sorted
3. Warn if A features are very large (>1Mb) and `-d` is used

## Common Patterns

**Basic coverage stats:**
```bash
bedtools coverage -a genes.bed -b reads.bed
```

**Coverage histogram:**
```bash
bedtools coverage -a genes.bed -b reads.bed -hist
```

**Mean coverage per feature:**
```bash
bedtools coverage -a genes.bed -b reads.bed -mean
```

**Memory-efficient sorted coverage:**
```bash
bedtools coverage -a sorted_genes.bed -b sorted_reads.bed -sorted
```

## Chaining Hints

- A features are often produced by `bedtools merge` or `bedtools slop`
- Output can be filtered with `awk` based on coverage fraction (column 7 in default output)
- Use `-hist` output with `grep "all"` to get genome-wide summary
