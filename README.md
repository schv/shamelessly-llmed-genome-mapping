# Genome Mapping Project

E. coli genome read mapping using suffix array-based algorithms.

## Requirements

- g++ with C++23 support
- ~200 MB disk space for data files

## Quick Start

### Build the Mapper

```bash
g++ -std=c++23 -O3 -o mapper mapper.cpp
```

### Run

```bash
# Map all reads (takes ~90 seconds)
./mapper

# Map first N reads (for quick testing)
./mapper -n 10000

# Custom parameters
./mapper -g data/genome.fna -r data/reads.fastq -n 100000 -s 20 -e 3
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-g <file>` | Reference genome (FASTA) | `data/GCF_000005845.2_ASM584v2_genomic.fna` |
| `-r <file>` | Reads file (FASTQ) | `data/ERR022075_1.fastq` |
| `-n <num>` | Max reads to process (-1 = all) | -1 |
| `-s <len>` | Seed length for mapping | 20 |
| `-e <num>` | Max edit distance allowed | 3 |
| `-h` | Show help | - |

## Data Files

Download and place in `data/` directory:

```bash
mkdir -p data && cd data

# Reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz

# Reads (FASTQ)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_1.fastq.gz
gunzip ERR022075_1.fastq.gz
```

## Library

The `lib/` folder contains reusable bioinformatics algorithms:

```cpp
#include "lib/bio.hpp"

// Suffix array
auto sa = bio::buildSuffixArray(text);
auto positions = bio::findAllOccurrences(text, sa, pattern);

// BWT
auto bwt = bio::computeBWT(text);
auto original = bio::inverseBWT(bwt);

// Edit distance
int dist = bio::editDistance<100>(s1, s2);

// K-mer analysis
auto [kmer, freq] = bio::findMostFrequentKmer(text, k);
```

## Output

The mapper outputs statistics including:
- Mapping rate (% reads mapped)
- Unique vs multi-mapped reads
- Average edit distance
- Genome coverage percentage
- Average sequencing depth

See `REPORT.md` for detailed results.
