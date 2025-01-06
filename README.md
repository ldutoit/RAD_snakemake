# README.md

This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline implementing the [Stacks](https://catchenlab.life.illinois.edu/stacks/manual/) 'denovo_map.pl' and 'refmap.pl' workflows for paired-end RAD-Sequencing data.

## Quick start

To run the pipeline on the example data (3 samples, denovo).

1.  Download this repository.
2.  From inside the repository, after checking dependencies:
   ```sh
  snakemake --cores all filtered.recode.vcf
   ```

## Dependencies

This pipeline was built with 

snakemake v 8.1.0
stacks v2.67
fastqc v0.12.1 
cutadapt v4.4
vcftools v0.1.5

for the reference-based version, the following dependencies are also required:

BWA v0.7.18
SAMTOOLS v1.9

## Tutorial

### Introduction

### Configuration

### Running the pipeline
