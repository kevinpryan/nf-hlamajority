# nf-hlamajority

![nf-hlamajority metro map](docs/images/test-bam-input-mosdepth-nf-hlamajority-with-subset-markdup.drawio.svg)

## Background

This pipeline is an implementation of a majority voting approach for the prediction of MHC Class I genotypes from DNA sequencing data. This approach was proposed by Claeys et al 2023 based on their benchmarking study. `nf-hlamajority` takes paired-end DNA-sequencing data and runs four tools:

- Optitype
- Polysolver
- Kourami
- HLA*LA

The MHC genotypes predicted by the highest number of tools is selected.

## Usage

Clone the repository

```bash
git clone https://github.com/kevinpryan/nf-hlamajority.git
```

### Build references

`--build_references` triggers a parallel workflow to build references, which is a prerequisite to running the pipeline. `--outdir` is the desired pipeline log directory. `--references_basedir` is the path to the directory where you wish to store your references. Default: `nf-hlamajority/references`.

```bash
nextflow run main.nf \
             --build_references \
             --outdir <PIPELINE_LOGS_OUTDIR> \
             --references_basedir <PATH_TO_REFERENCES_DIR> \
             -profile <singularity/docker/awsbatch>
```

This workflow carries out the following steps:

- Build BWAkit references
- Download the Kourami reference
- Build Kourami database
- Download HLA*LA reference
- Compute HLA*LA graph index structure
- Index reference FASTAs from BWAkit, Kourami, HLA*LA

A local test of the reference building workflow on a SLURM HPC using Singularity took 2 hours 43 minutes to run, and required a maximum of 33.4 GB of RAM. This reference is static and can be reused across genotyping runs.

Your references directory should have the following structure:

```bash
references
├── bwakit
│   ├── hs38DH.fa
│   ├── hs38DH.fa.alt
│   ├── hs38DH.fa.amb
│   ├── hs38DH.fa.ann
│   ├── hs38DH.fa.bwt
│   ├── hs38DH.fa.pac
│   └── hs38DH.fa.sa
├── hla-la
│   └── PRG_MHC_GRCh38_withIMGT
├── kourami
│   ├── build.xml
│   ├── custom_db
│   ├── LICENSE
│   ├── pom.xml
│   ├── preprocessing.md
│   ├── README.md
│   ├── resources
│   ├── scripts
│   ├── src
│   └── target
└── source
    └── IMGTHLA
```

When you have installed and built the required references, you are ready to run the pipeline.

### Running nf-hlamajority

The pipeline accepts the following input file types:

- FASTQ
- Aligned BAM
- CRAM

It is designed for paired-end DNA sequencing data, but will also accept single-end data. However, HLA genotyping using single-end data is less reliable than paired-end and is not recommended.

#### Test data

To ensure the pipeline is working as expected, the test profile should be run first.

```bash
nextflow run main.nf \
       --outdir <OUTDIR> \
       -profile <test,singularity/docker/awsbatch>
```

The test dataset is a CRAM file (316 MB) provided through the [HLA*LA GitHub repository](https://github.com/DiltheyLab/HLA-LA/tree/master) (1000 Genomes sample NA12878).

The expected outputs of each tool from the test dataset can be found at `assets/test-outputs/test-outputs-1000genomes/NA12878/`

#### Running on full datasets
       
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

When using aligned data, you can provide a samplesheet containing both BAM and CRAM files. They do not need to be sorted or indexed; coordinate sorting is performed internally where required.

You must pass the `--aligned` flag when using BAM or CRAM files as input.

When using CRAM files, you must pass the reference fasta used to generate the CRAM file via the `--cram_fasta` parameter. The pipeline only supports one `--cram_fasta` per run.

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

By default, the pipeline uses the majority voting method proposed by Claeys et al, whereby each tool gets one vote and the genotype with the most votes is assigned. In the case of a tie, the genotype of the best-performing tool in the benchmark is assigned (`--voting_method majority`).

An alternative is to carry out a weighted vote (`--voting_method weighted`). By default, the pipeline uses the accuracy scores for each tool in the Claeys et al benchmark (for each HLA gene) as the weight (`assets/benchmarking_results_claeys_cleaned.csv`). Weighted voting prioritises tools that demonstrated higher per-gene accuracy in the Claeys et al. benchmark, allowing higher-confidence calls to dominate in cases of disagreement.

The user can specify their own weights by providing their own csv file to `--weights` in the following format:

```bash
tool,A,B,C
hlala,0.899,0.972,0.962
kourami,0.834,0.761,0.796
optitype,0.98,0.976,0.984
polysolver,0.949,0.918,0.98
```
so for example

*for FASTQ output using weighted voting with the default weights:*

```bash
nextflow run main.nf \
       --samplesheet <SAMPLESHEET> \
       --outdir <OUTDIR> \
       --voting_method weighted \
       -profile <singularity/cluster/.../institute>
```

Whatever method is used, the following cross-sample output files are expected:

```bash
├── nf_hlamajority_all_calls_sorted.tsv
├── nf_hlamajority_depth_sorted.tsv
├── nf_hlamajority_stats_combined_sorted.tsv
└── nf_hlamajority_votes_combined_sorted.tsv
```

Example outputs:

*nf_hlamajority_all_calls_sorted.tsv*

```bash
sample	tool	A1	A2	B1	B2	C1	C2
NA12878	hlala	11:01	01:01	08:01	56:01	01:02	07:01
NA12878	kourami	01:01	11:01	56:01	08:01	01:02	07:01
NA12878	optitype	01:01	11:01	08:01	56:01	01:02	07:01
NA12878	polysolver	01:01	11:01	08:01	56:01	01:02	07:01
```

*nf_hlamajority_depth_sorted.tsv*

```bash
sample	gene	mean_depth_hla_exons_2_3_gene	mean_depth_hla_exons_2_3_classI
NA12878	HLA-A	15.9	19.79
NA12878	HLA-B	19.74	19.79
NA12878	HLA-C	23.73	19.79
```

*nf_hlamajority_stats_combined_sorted.tsv*

```bash
sample	gene	allele1	allele2	support	matching_tools	method	weight_winner	total_weight	n_tools_support	n_tools_called	mean_depth_hla_exons_2_3_gene	mean_depth_hla_exons_2_3_classI
NA12878	HLA-A	01:01	11:01	1	hlala,kourami,optitype,polysolver	majority_vote	4	4	4	4	15.9	19.79
NA12878	HLA-B	08:01	56:01	1	hlala,kourami,optitype,polysolver	majority_vote	4	4	4	4	19.74	19.79
NA12878	HLA-C	01:02	07:01	1	hlala,kourami,optitype,polysolver	majority_vote	4	4	4	4	23.73	19.79
```

*nf_hlamajority_votes_combined_sorted.tsv*

```bash
sample	gene	allele1	allele2	matching_tools	method	support	mean_depth_hla_exons_2_3_gene
NA12878	HLA-A	01:01	11:01	hlala,kourami,optitype,polysolver	majority_vote	1	15.9
NA12878	HLA-B	08:01	56:01	hlala,kourami,optitype,polysolver	majority_vote	1	19.74
NA12878	HLA-C	01:02	07:01	hlala,kourami,optitype,polysolver	majority_vote	1	23.73
```

## Dependencies

The pipeline requires:
- Nextflow (DSL2)
- Singularity/Apptainer or Docker
- Java (compatible with your Nextflow version)


## Integration with Landscape of Effective Neoantigens Software (LENS)

Instructions for running nf-hlamajority as part of the LENS neoantigen prediction pipeline can be found at `docs/lens-hlamajority/`.

 ## Important Licensing Information

This pipeline is licensed under the MIT license. However, it is designed to be run with a Docker/Singularity container that includes the POLYSOLVER software. The use of this specific POLYSOLVER container is subject to its own licensing terms, which are more restrictive.

As per the POLYSOLVER license:

- Academic and Non-Profit Use: The container is available for academic and non-profit users.
- Commercial Use: Commercial use of this pipeline with the provided POLYSOLVER container is prohibited unless you have first obtained a commercial license for Novoalign and any other applicable software within the container.

By using this pipeline, you are agreeing to comply with the licensing terms of all its software components. It is the user's responsibility to ensure they are licensed appropriately. For commercial use, please contact Novocraft Technologies to acquire a license for Novoalign.

## References

Claeys, A., Merseburger, P., Staut, J., Marchal, K., & van den Eynden, J. (2023). Benchmark of tools for in silico prediction of MHC class I and class II genotypes from NGS data. BMC Genomics, 24(1), 1–14. https://doi.org/10.1186/s12864-023-09351-z
