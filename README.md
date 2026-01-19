# nf-hlamajority

![nf-hlamajority metro map](docs/images/nf-hlamajority-with-subset-markdup.drawio.svg)

## Background

This pipeline is an implementation of a majority voting approach for the prediction of MHC Class I genotypes from DNA sequencing data. This approach was proposed by Claeys et al 2023 based on their benchmarking study. `nf-hlamajority` takes paired-end DNA-sequencing data and runs four tools:

- Optitype
- Polysolver
- Kourami
- HLA-LA

The MHC genotypes predicted by the highest number of tools is chosen.

This pipeline contains a step where the aligned BAM is subset to chromosome 6, chromosome 6 alt contigs and HLA contigs. Claeys et al did not subset BAMs in this way in their benchmark. In addition, only a version of this Nextflow pipeline that does not subset the BAMs has been tested on a portion of the data used in the Claeys et al benchmark. Future work will include repeating these tests for the version of the pipeline that subsets the aligned BAMs.

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

The script `install_references.sh` must be run before running the pipeline for the first time. The script should be run from **within** the `bin` directory.

The script takes an argument which is a directory (which must exist) where you can download the singularity images to be used to build/download the references.

```bash
cd bin
bash install_references.sh /path/to/singularity/cache/dir/
```

This will take several hours to run, and can require up to ~40Gb memory.


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

Claeys, A., Merseburger, P., Staut, J., Marchal, K., & van den Eynden, J. (2023). Benchmark of tools for in silico prediction of MHC class I and class II genotypes from NGS data. BMC Genomics, 24(1), 1–14. https://doi.org/10.1186/s12864-023-09351-z
