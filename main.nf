#!/usr/bin/env nextflow

include { HLATYPING } from "./workflows/hlatyping"

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
        params.reference_basename,
        params.reference_dir,
        //params.hlatypes,
        params.chr,
        params.hla_la_graph,
        params.kourami_ref,
        params.kourami_database,
        params.trimmer,
        params.save_trimmed_fail,
        params.save_merged//,
        //params.adapter_fasta,
        //params.subset_regions
    )
}

// if you want to add back in seq_type option to add rna-seq:
/*
Channel.fromPath(params.samplesheet, checkIfExists: true)
    | splitCsv( header:true )
    | map { row ->
        meta = row.subMap('sample', 'seq_type')
        [meta, [
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)]]
    }
    | set { ch_fastq }
*/
