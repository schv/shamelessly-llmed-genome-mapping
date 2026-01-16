# Genome Mapping Report

## Overview
This report presents the results of mapping Illumina sequencing reads to the *Escherichia coli* str. K-12 substr. MG1655 reference genome.

## Data

### Reference Genome
- **Source**: NCBI (GCF_000005845.2_ASM584v2)
- **Organism**: *Escherichia coli* str. K-12 substr. MG1655
- **Genome Size**: 4,641,652 bp

### Sequencing Reads
- **Source**: EBI SRA (ERR022075)
- **Format**: FASTQ (Illumina)
- **Read Length**: 100 bp
- **Total Reads**: 22,720,100

## Algorithms Used

1. **Suffix Array** - construction for efficient pattern matching on the reference genome
2. **Seed-and-Extend** - 20-mer seeds extracted from reads and searched in the suffix array
3. **Band-Limited Edit Distance** - Verification of candidate alignments allowing up to 3 mismatches/indels

### Implementation Details
- Suffix array built once for the reference genome (~1 second)
- Each read is first tested for exact match (fast path)
- If no exact match, seeds are extracted and candidate positions verified with edit distance
- Multi-mapping detection for reads matching multiple genome locations

## Results

### Mapping Statistics

| Metric | Count | Percentage |
|--------|-------|------------|
| Total Reads Processed | 22,720,100 | 100.00% |
| Mapped Reads | 11,219,159 | 49.38% |
| Unmapped Reads | 11,500,941 | 50.62% |
| Uniquely Mapped | 10,897,207 | 47.96% |
| Multi-Mapped | 321,952 | 1.42% |

### Alignment Quality

| Metric | Value |
|--------|-------|
| Average Edit Distance | 0.16 |
| Alignments with 0 errors | ~84% |
| Max Allowed Errors | 3 |

### Genome Coverage (from uniquely mapped reads)

| Metric | Value |
|--------|-------|
| Covered Bases | 4,571,846 bp |
| Coverage Percentage | 98.50% |
| Average Depth | 234.77x |
| Uncovered Bases | 69,806 bp (1.50%) |
