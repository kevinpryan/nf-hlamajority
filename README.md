# nf-hlamajority

![nf-hlamajority metro map](docs/images/nf-hlamajority-without-subset.svg)

## Background

This pipeline is an implementation of a majority voting approach for the predition of MHC Class I genotypes from DNA sequencing data. This approach was proposed by Claeys et al 2023 based on their benchmarking study. `nf-hlamajority` takes paired-end DNA-sequencing data and runs four tools:

- Optitype
- Polysolver
- Kourami
- HLA-LA

The MHC genotypes predicted by the highest number of tools is chosen.

## Usage

Clone the repository

```bash
git clone https://github.com/kevinpryan/nf-hlamajority.git
```

Prepare a samplesheet with your input data that looks as follows:

**samplesheet.csv**

```csv
sample,fastq_1,fastq_2
SAMPLE1,SAMPLE1_S1_L002_R1_001.fastq.gz,SAMPLE1_S1_L002_R2_001.fastq.gz
SAMPLE2,SAMPLE2_S1_L003_R1_001.fastq.gz,SAMPLE2_S1_L003_R2_001.fastq.gz
SAMPLE3,SAMPLE3_S1_L004_R1_001.fastq.gz,SAMPLE3_S1_L004_R2_001.fastq.gz
```

Each row represents a pair of fastq files. Currently the pipeline only supports paired-end fastqs.

When you have installed and built the required references, run the pipeline with:

```bash
nextflow run main.nf \
       --samplesheet <SAMPLESHEET> \
       --outdir <OUTDIR> \
       -profile <singularity/cluster/.../institute>
```

## Dependencies

The pipeline requires the following:

- nextflow
- singularity

## References
