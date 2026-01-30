---
name: bedtools-intersect
description: Find overlapping intervals between two genomic feature files using bedtools intersect.
---

# bedtools intersect

## Purpose

Find overlapping intervals between two feature files (BED/GFF/VCF/BAM).

## Key Arguments

| Flag | Description |
|------|-------------|
| `-a` | Input file A (query) |
| `-b` | Input file(s) B (database). Multiple files and wildcards allowed |
| `-wa` | Write original A entry for each overlap |
| `-wb` | Write original B entry for each overlap |
| `-wo` | Write A and B entries plus overlap base pairs |
| `-wao` | Like `-wo` but also reports A entries with no overlap (overlap=0) |
| `-loj` | Left outer join: report all A entries, NULL B if no overlap |
| `-v` | Report A entries with NO overlap in B (like grep -v) |
| `-u` | Report each A entry once if any overlap found |
| `-c` | For each A entry, report count of overlaps with B |
| `-f` | Minimum overlap as fraction of A (default: 1bp) |
| `-F` | Minimum overlap as fraction of B |
| `-r` | Require reciprocal overlap (both `-f` and `-F` must be met) |
| `-e` | Require that the minimum fraction be satisfied for A OR B (instead of AND) |
| `-s` | Require same strandedness |
| `-S` | Require opposite strandedness |
| `-sorted` | Use sweep algorithm (requires sorted input, uses less memory) |
| `-header` | Print header from A file before results |

## Assumptions

- Both `-a` and `-b` files are tab-delimited with at least 3 columns (chrom, start, end)
- Chromosome names are consistent across files (e.g., both use `chr1` or both use `1`)
- BED coordinates are 0-based, half-open
- If using `-sorted`, both files must be sorted by chromosome then start position
- If using `-s` or `-S`, files must have a strand column (column 6 in BED)

## Common Issues

- **Empty output**: Most often caused by chromosome naming mismatch (`chr1` vs `1`). Check with: `cut -f1 file | sort -u | head`
- **Using `-sorted` on unsorted input**: bedtools will silently produce wrong results or error. Always verify sort order first
- **Off-by-one confusion**: BED is 0-based half-open. `chr1 10 20` covers bases 10-19. Two intervals like `chr1 10 20` and `chr1 20 30` do NOT overlap
- **Unexpected duplicates**: If A has overlapping entries and B has overlapping entries, you may get more output rows than expected. Use `-u` to deduplicate A entries
- **Forgetting `-wa`/`-wb`**: Default output is the intersection interval itself, not the original entries. Use `-wa` and/or `-wb` to preserve original entries

## Validation

Before running:
1. Confirm both A and B files exist and have >= 3 tab-delimited columns
2. If `-sorted` is requested, verify both files are sorted
3. If `-s`/`-S` is used, verify strand column exists (>= 6 columns in BED)
4. Check chromosome naming consistency between A and B

## Common Patterns

**Basic overlap — which A features overlap B:**
```bash
bedtools intersect -a peaks.bed -b genes.bed -wa
```

**Exclude blacklisted regions:**
```bash
bedtools intersect -a peaks.bed -b blacklist.bed -v
```

**Reciprocal 50% overlap:**
```bash
bedtools intersect -a a.bed -b b.bed -f 0.5 -r
```

**Count overlaps per feature:**
```bash
bedtools intersect -a genes.bed -b reads.bed -c
```

**Report both entries and overlap size:**
```bash
bedtools intersect -a a.bed -b b.bed -wo
```

**Memory-efficient sorted intersection:**
```bash
bedtools intersect -a sorted_a.bed -b sorted_b.bed -sorted
```

## Chaining Hints

- Output is often piped to `bedtools merge` to collapse overlapping results
- Commonly preceded by `bedtools sort` when using `-sorted`
- Use `-v` output as input to `bedtools closest` to find nearest non-overlapping feature
- Pipe to `bedtools summary` for QC on the result
