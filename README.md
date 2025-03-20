# SARS-CoV-2 Variant Analysis Pipeline

A complete bioinformatics pipeline for analyzing SARS-CoV-2 genomic data, with a focus on variant calling and annotation. This pipeline can be run using either Bash scripts or Snakemake for workflow management.

## Overview

This pipeline performs the following analysis steps:

1. Downloading FASTQ data from SRA using SRA-toolkit
2. Quality control of raw sequences using FastQC
3. Adapter and quality trimming with Trimmomatic
4. Quality control of trimmed sequences
5. Generating MultiQC reports for quality assessment
6. Alignment to the SARS-CoV-2 reference genome using BWA-MEM
7. Converting, sorting, and indexing BAM files with SAMtools
8. Variant calling using FreeBayes
9. Variant filtering with GATK
10. Variant annotation with SnpEff
11. Generating readable variant reports

## Repository Structure

```
sars-cov2-variant-analysis/
├── README.md                      # Project documentation
├── environment.yml                # Conda environment file
├── data/                          # Input data directory
│   ├── SARS-CoV-2-reference.fasta # Reference genome
│   └── adapters.fa                # Adapter sequences
├── bash/                          # Bash pipeline implementation
│   ├── pipeline.sh                # Main pipeline script
│   └── run_analysis.sh            # Wrapper script with CLI options
├── snakemake/                     # Snakemake pipeline implementation
│   ├── Snakefile                  # Workflow definition
│   ├── config.yaml                # Configuration parameters
│   └── run_snakemake.sh           # Snakemake execution script
├── docker/                        # Docker configuration
│   ├── Dockerfile                 # Container definition
│   ├── build_docker.sh            # Docker build script
│   └── run_docker.sh              # Docker run script
├── results/                       # Output directory (created during execution)
└── docs/                          # Documentation and project writeup
    └── images/                    # Images for documentation
```

## Requirements

All required software is specified in the `environment.yml` file and can be installed with Conda/Mamba:

```bash
conda env create -f environment.yml
```

Or you can use the provided Docker container which includes all dependencies:

```bash
cd docker
./build_docker.sh
./run_docker.sh
```

## Usage

### Bash Pipeline

```bash
cd bash
./run_analysis.sh
```

For more options:

```bash
./run_analysis.sh --help
```

### Snakemake Pipeline

```bash
cd snakemake
./run_snakemake.sh
```

For more options:

```bash
./run_snakemake.sh --help
```

## Results

The analysis pipeline creates the following key output files:

- Quality control reports in `results/qc_results/`
- Trimmed reads in `results/trimmed_data/`
- Alignment files (SAM/BAM) in `results/`
- Variant calls (VCF) in `results/`
- Annotated variants in `results/annotated_variants.vcf`
- Readable variant summary in `results/gene_annotations.txt`

## Citation

If you use this pipeline in your research, please cite:

```
Lundquist, K. (2025). SARS-CoV-2 Variant Analysis Pipeline. GitHub. 
https://github.com/yourusername/sars-cov2-variant-analysis
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.