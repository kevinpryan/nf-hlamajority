#!/usr/bin/env nextflow

include { HLATYPING } from "./workflows/hlatyping"
include { samtools_sort } from "./modules/local/samtools_sort"
include { samtools_index } from "./modules/local/samtools_index"
include { subsetBam2 } from "./modules/local/subsetBam"
include { bam2fastq } from "./modules/local/bam2fastq"

workflow {
    Channel.fromPath(params.samplesheet, checkIfExists: true)
    | splitCsv( header:true )
    | map { row ->
        meta = row.subMap('sample')
        [meta, [
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)]]
    }
    | set { ch_fastq }
// example ch_fastq: [[sample:3532, seq_type:dna], [/data4/kryan/misc/useful/nextflow/nf-hlatyping/testdir/gen_testdata/3532_subset_10000.1.fq.gz, /data4/kryan/misc/useful/nextflow/nf-hlatyping/testdir/gen_testdata/3532_subset_10000.2.fq.gz]]
    HLATYPING(
        ch_fastq,
        params.reference_dir,
        params.hla_la_graph,
        params.kourami_ref,
        params.kourami_database,
        params.trimmer,
        params.adapter_fasta,
        params.save_trimmed_fail,
        params.save_merged//,
    )
    
}
