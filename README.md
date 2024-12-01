# Bifidotyper

Bifidotyper is a bioinformatics tool designed to take you from raw FastQ files to a complete and reproducible analysis of [Bifidobacteria](https://en.wikipedia.org/wiki/Bifidobacterium) strains in your samples. It makes extensive use of [Sylph](https://www.nature.com/articles/s41587-024-02412-y) for rapid, k-mer-based read alignments. It also detects the presence of genes necessary for the metabolism of human milk oligosaccharides ([HMOs](https://en.wikipedia.org/wiki/Human_milk_oligosaccharide)) based on alignments to the Bifidobacterium longum genome ([CP001095.1](https://www.ncbi.nlm.nih.gov/nuccore/CP001095.1/)) using gene annotations from [Brodin et al.](https://doi.org/10.1016/j.cell.2021.05.030).

Bifidotyper was developed as part of a PhD rotation in the [Olm Lab](https://www.colorado.edu/lab/olm/).

## Features

- Sketch genomes and reads using Sylph.
- Query and profile genomes.
- Detect HMO genes in sequencing data.
- Generate informative plots for taxonomic abundance and containment indices.

## Installation

Clone the repository and install the package:
```bash
git clone https://github.com/Bennibraun/bifidotyper.git
cd bifidotyper
```

Create a Conda environment for required software:
```bash
conda env create --name bifidotyper --file=env.yaml
```

Then you can simply run:
```bash
pip install -e .
```


## Usage

### Command Line Interface

For single-end reads:
```bash
bifidotyper -se <single-end FASTQ files> [options]
```

Or paired-end reads:
```bash
bifidotyper -pe <paired-end FASTQ files> [options]
```

### Options

- `-se, --single-end`: Single-end FASTQ files.
- `-pe, --paired-end`: Paired-end FASTQ files (R1 and R2 files, supports wildcards).
- `--r1-suffix`: Suffix for R1 files (only for paired-end mode).
- `--r2-suffix`: Suffix for R2 files (only for paired-end mode).
- `-t, --threads`: Number of threads to use for parallel processing (default: 1).

### Example
```bash
bifidotyper -pe data/*.fastq.gz --r1-suffix _R1 --r2-suffix _R2 -t 4
```

## Output

The tool generates several output files and directories:

- `sylph_genome_sketches/`: Contains genome sketches.
- `sylph_fastq_sketches/`: Contains FASTQ sketches.
- `sylph_genome_queries/`: Contains genome query and profile results.
- `hmo_quantification/`: Contains HMO gene quantification results.
- `plots/`: Contains generated plots.

## Logging

A log file `bifidotyper.log` is created in the working directory, capturing detailed information about the processing steps. Specific logs for each step are including in their respective output directories.

---

# Files
- humann2_HMO_annotation.csv was retrieved from [Brodin et al. 2021](https://data.mendeley.com/datasets/gc4d9h4x67/2)
- All B. longum annotations are from the NCBI record for [CP001095.1](https://www.ncbi.nlm.nih.gov/nuccore/CP001095.1/)
- For simplicity, this release currently includes all required Bifidobacterial genomes. In the future, a module will be included to download an updated set of reference genomes as part of the program's execution.

---

## License

This project is licensed under the MIT License.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.
