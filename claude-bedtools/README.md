# claude-bedtools

Claude skills for genomic interval analysis using [bedtools](https://bedtools.readthedocs.io/).

Describe genomic analyses in natural language and Claude will run the right bedtools commands with proper arguments, input validation, and quality control.

## Skills

| Skill | Description |
|-------|-------------|
| `bedtools` | Parent dispatcher — detects genomic analysis tasks, validates inputs, orchestrates pipelines |
| `bedtools-intersect` | Find overlapping intervals between files |
| `bedtools-merge` | Combine overlapping/nearby intervals |
| `bedtools-sort` | Sort intervals by position or other criteria |
| `bedtools-closest` | Find nearest feature in another file |
| `bedtools-coverage` | Compute depth and breadth of coverage |
| `bedtools-genomecov` | Compute genome-wide coverage |
| `bedtools-subtract` | Remove overlapping portions of intervals |
| `bedtools-complement` | Find uncovered regions of the genome |
| `bedtools-slop` | Extend intervals by a specified amount |
| `bedtools-getfasta` | Extract DNA sequences for intervals |
| `bedtools-summary` | QC and summary statistics |

## Requirements

- [bedtools](https://bedtools.readthedocs.io/) installed and in PATH
- Claude Code with skills support

## Installation

Copy or symlink the `skills/` directory into your project or `~/.claude/` directory.

## Usage

Simply describe what you want to do:

- "Find which peaks overlap with known genes"
- "Merge overlapping intervals and count how many were combined"
- "What is the closest gene to each of my variants?"
- "Compute genome-wide coverage from my BAM file"
- "Extract sequences for these regions"

The skill handles input validation, chromosome naming checks, sort order verification, and error recovery automatically.
