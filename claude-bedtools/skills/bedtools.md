---
name: bedtools
description: Genomic interval analysis using bedtools. Use when users mention bedtools, genomic overlaps/intersections, merging regions, nearest features, coverage, complement, subtraction, or working with BED/BAM/VCF/GFF files. Also triggers on domain tasks like "find peaks overlapping genes" or "compute genome coverage".
---

# Bedtools — Genomic Interval Analysis

You are guiding Claude to perform genomic interval analysis using bedtools.

## When This Skill Applies

Invoke this skill when the user:
- Mentions bedtools by name
- Describes genomic interval operations: overlaps, intersections, merging regions, nearest features, coverage, complement, subtraction
- References file types: BED, BAM, VCF, GFF/GTF, BEDPE
- Describes domain tasks: "find peaks overlapping genes", "compute genome coverage", "extract sequences from regions", "find nearest gene", "remove blacklisted regions"

## Pre-flight Checks (ALWAYS run these first)

Before executing any bedtools command, perform these checks:

### 1. Verify bedtools installation
Run: `which bedtools && bedtools --version`
If not found, inform the user that bedtools must be installed and suggest: `conda install -c bioconda bedtools` or `brew install bedtools`.

### 2. Validate input files exist
For each input file the user references, run: `test -f <file> && echo "OK" || echo "NOT FOUND: <file>"`

### 3. Sniff file format
For each input file, run:
```bash
head -n 5 <file> | awk -F'\t' '{print NF}' | sort -u
```
Check:
- Files are tab-delimited
- BED files have at least 3 columns
- Look for header lines (`track`, `browser`, lines starting with `#`) that may need `-header` flag
- Warn if column count is inconsistent across lines

### 4. Check chromosome naming consistency
When multiple input files are used, compare chromosome names:
```bash
cut -f1 <file1> | sort -u | head -20 > /tmp/chroms1.txt
cut -f1 <file2> | sort -u | head -20 > /tmp/chroms2.txt
comm -3 /tmp/chroms1.txt /tmp/chroms2.txt
```
Warn if one file uses `chr1` and the other uses `1`.

### 5. Check sort order (when required)
For commands that require sorted input (`merge`, `closest`, `complement`):
```bash
bedtools sort -i <file> | diff - <file> | head -5
```
If output is non-empty, warn the user and offer to sort first.

## Dispatch to Subcommand Skills

Based on the user's intent, invoke the appropriate subcommand skill:

| User intent | Skill |
|---|---|
| Find overlapping intervals | `bedtools-intersect` |
| Combine/merge overlapping intervals | `bedtools-merge` |
| Sort intervals | `bedtools-sort` |
| Find nearest feature | `bedtools-closest` |
| Compute depth/breadth of coverage | `bedtools-coverage` |
| Compute genome-wide coverage | `bedtools-genomecov` |
| Remove overlapping portions | `bedtools-subtract` |
| Find uncovered regions | `bedtools-complement` |
| Extend/pad intervals | `bedtools-slop` |
| Extract DNA sequences | `bedtools-getfasta` |
| QC / summary statistics | `bedtools-summary` |

## Chaining and Pipelines

When the user's request implies multiple steps:
1. Plan the full pipeline before executing
2. Present the plan to the user
3. Execute step by step
4. Use pipes by default; write to named files for multi-step pipelines where intermediate inspection is useful
5. Report row counts at each step

Example: "find merged peaks that overlap promoters"
1. `bedtools sort -i peaks.bed`
2. `bedtools merge -i -` (piped)
3. `bedtools intersect -a - -b promoters.bed`

## Output Handling

After running commands:
1. Report row count of output: `wc -l < output.bed`
2. Flag potential issues:
   - Empty output (0 lines) — likely a chromosome naming mismatch or format problem
   - Unexpectedly large output — may indicate overly broad matching
3. Offer to run `bedtools summary` for QC and summary analysis
4. For complex pipelines, always suggest `bedtools summary` on the final output

## Error Recovery

If a bedtools command fails:
1. Read stderr carefully
2. Diagnose the issue (missing file, wrong format, unsorted input, missing genome file)
3. Suggest a fix
4. If the fix is a pre-processing step (e.g., sorting), offer to add it to the pipeline
