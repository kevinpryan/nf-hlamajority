# nf-hlamajority

![nf-hlamajority metro map](docs/images/bam-input-mosdepth-nf-hlamajority-with-subset-markdup.drawio.svg)

## Background

This pipeline is an implementation of a majority voting approach for the prediction of MHC Class I genotypes from DNA sequencing data. This approach was proposed by Claeys et al 2023 based on their benchmarking study. `nf-hlamajority` takes paired-end DNA-sequencing data and runs four tools:

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

The pipeline accepts the following input file types:

- FASTQ
- Aligned BAM
- CRAM

It is designed for paired-end sequencing data.

Prepare a samplesheet with your input data that looks as follows:

**samplesheet.csv**

```csv
sample,fastq_1,fastq_2
SAMPLE1,SAMPLE1_S1_L002_R1_001.fastq.gz,SAMPLE1_S1_L002_R2_001.fastq.gz
SAMPLE2,SAMPLE2_S1_L003_R1_001.fastq.gz,SAMPLE2_S1_L003_R2_001.fastq.gz
SAMPLE3,SAMPLE3_S1_L004_R1_001.fastq.gz,SAMPLE3_S1_L004_R2_001.fastq.gz
```

or if you are using an aligned data type (BAM, CRAM), prepare the samplesheet as follows:

```csv
sample,aln
SAMPLE1,SAMPLE1.bam
SAMPLE2,SAMPLE2.cram
```

When using aligned data, you can provide a samplesheet containing both BAM and CRAM files. They do not need to be sorted or indexed.

You must pass the `--aligned` flag when using BAM or CRAM files as input.

When using CRAM files, you must pass the reference fasta used to generate the CRAM file via the `--cram_fasta` parameter. The pipeline only supports one `--cram_fasta` per run.

The script `install_references.sh` must be run before running the pipeline for the first time. The script should be run from **within** the `bin` directory.

The script takes an argument which is a directory (which must exist) where you can download the singularity images to be used to build/download the references.

```bash
cd bin
bash install-references.sh --engine <docker,singularity> --cache /path/to/singularity/cache/dir/
```

`--cache` is only required when using singularity.

A local test of `install-references.sh` on a SLURM HPC using singularity took 1 hour 40 minutes to run, and required approximately 33.4 GB of RAM. This reference is static and can be reused across genotyping runs.

When you have installed and built the required references, run the pipeline with:

*for FASTQ input*

```bash
nextflow run main.nf \
       --samplesheet <SAMPLESHEET> \
       --outdir <OUTDIR> \
       -profile <singularity/cluster/.../institute>
```

*for BAM input*

```bash
nextflow run main.nf \
       --samplesheet <SAMPLESHEET> \
       --outdir <OUTDIR> \
       --aligned \
       -profile <singularity/cluster/.../institute>
```

*for aligned input including at least one CRAM*

```bash
nextflow run main.nf \
       --samplesheet <SAMPLESHEET> \
       --outdir <OUTDIR> \
       --aligned \
       --cram_fasta \
       -profile <singularity/cluster/.../institute>
```

## Dependencies

The pipeline requires the following:

- nextflow
- singularity or docker

## References

Claeys, A., Merseburger, P., Staut, J., Marchal, K., & van den Eynden, J. (2023). Benchmark of tools for in silico prediction of MHC class I and class II genotypes from NGS data. BMC Genomics, 24(1), 1–14. https://doi.org/10.1186/s12864-023-09351-z
