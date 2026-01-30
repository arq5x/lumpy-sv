# Bedtools Skill Test Scenarios

Test these natural language prompts to verify the skills work correctly.

## Basic Dispatch Tests

1. **Explicit bedtools mention**: "Run bedtools intersect on a.bed and b.bed"
   - Expected: triggers bedtools skill, runs pre-flight, dispatches to intersect

2. **Natural language overlap**: "Find which of my peaks overlap with known genes"
   - Expected: triggers bedtools skill, dispatches to intersect

3. **Merge request**: "Combine overlapping intervals in my BED file"
   - Expected: triggers bedtools skill, dispatches to sort then merge

4. **Nearest feature**: "What is the closest gene to each of my variants?"
   - Expected: triggers bedtools skill, dispatches to closest

## Validation Tests

5. **Missing file**: "Intersect nonexistent.bed with b.bed"
   - Expected: pre-flight catches missing file, reports error

6. **Chromosome mismatch**: Create a file with `1` instead of `chr1`, intersect with `chr1` file
   - Expected: pre-flight warns about naming mismatch

7. **Unsorted input for merge**: "Merge overlapping intervals in unsorted.bed"
   - Expected: detects unsorted input, offers to sort first

## Pipeline Tests

8. **Multi-step**: "Find merged peaks that overlap promoters"
   - Expected: plans sort → merge → intersect pipeline

9. **With QC**: "Intersect peaks with genes and give me a summary"
   - Expected: runs intersect then bedtools summary

## Edge Cases

10. **Empty result**: Intersect files with no overlaps
    - Expected: warns about empty output, suggests checking chromosome names
