---
name: bedtools-getfasta
description: Extract DNA sequences from a FASTA file for given intervals using bedtools getfasta.
---

# bedtools getfasta

## Purpose

Extract DNA sequences from a FASTA file based on feature coordinates.

## Key Arguments

| Flag | Description |
|------|-------------|
| `-fi` | Input FASTA file |
| `-bed` | BED/GFF/VCF file of regions to extract |
| `-fo` | Output file (default: stdout) |
| `-name` | Use name field and coordinates for FASTA header |
| `-nameOnly` | Use only the name field for FASTA header |
| `-tab` | Output in tab-delimited format instead of FASTA |
| `-bedOut` | Output in BED format with sequence as extra column |
| `-s` | Reverse complement if feature is on antisense strand |
| `-split` | For BED12, extract and concatenate block sequences (e.g., exons) |
| `-fullHeader` | Use full FASTA header (not just first word) |
| `-rna` | Treat as RNA (affects reverse complement) |

## Assumptions

- FASTA file must be indexed (`.fai` file must exist). Create with: `samtools faidx ref.fa`
- BED coordinates must not exceed chromosome lengths in the FASTA
- Chromosome names in BED must match FASTA headers exactly
- BED is 0-based, half-open

## Common Issues

- **Missing FASTA index**: `getfasta` requires a `.fai` index. Run `samtools faidx ref.fa` first
- **Chromosome name mismatch**: FASTA may use `>chr1` while BED uses `1` or vice versa
- **Coordinates out of bounds**: BED end position exceeds chromosome length. Check with: `awk -v OFS='\t' 'NR==FNR{len[$1]=$2;next} $3>len[$1]' ref.fa.fai regions.bed`
- **Missing `-s` flag**: If you need strand-aware sequences (e.g., for gene sequences), you must use `-s`. Without it, all sequences are from the forward strand
- **BED12 exon extraction**: Use `-split` to get spliced transcript sequences from BED12 entries

## Validation

Before running:
1. Confirm FASTA file exists and has `.fai` index
2. Confirm BED file exists
3. Check chromosome naming consistency between BED and FASTA
4. Optionally check that no BED coordinates exceed chromosome lengths

## Common Patterns

**Extract sequences for regions:**
```bash
bedtools getfasta -fi ref.fa -bed regions.bed
```

**Strand-aware extraction:**
```bash
bedtools getfasta -fi ref.fa -bed genes.bed -s -name
```

**Tab-delimited output:**
```bash
bedtools getfasta -fi ref.fa -bed regions.bed -tab
```

**Extract spliced exon sequences from BED12:**
```bash
bedtools getfasta -fi ref.fa -bed transcripts.bed12 -split -s -name
```

## Chaining Hints

- BED input often comes from `bedtools slop` (to get flanking sequences) or `bedtools intersect`
- Output FASTA can be piped to sequence analysis tools (MEME, BLAST, etc.)
- Use `-tab` or `-bedOut` for downstream parsing
