#!/usr/bin/env nextflow

include { HLATYPING } from "./workflows/hlatyping"

workflow {
    if (params.bam = true) {
    // parse BAM samplesheet
    Channel.fromPath(params.samplesheet, checkIfExists: true)
    | splitCsv( header:true )
    | map { row ->
        meta = row.subMap('sample')
        [meta, [
            file(row.bam, checkIfExists: true),
    }
    | set { ch_bam }
    // subset BAM to regions of interest
    samtools_sort(ch_bam)
    samtools_index(samtools_sort.out.sortedbam)
    subsetBam2(samtools_index.out.bam_indexed)
    bam2fastq(subsetBam2.out.subsetbam)
    // run hla typing pipeline
    HLATYPING(
        bam2fastq.out.convertedfastqs,
        params.reference_dir,
        params.hla_la_graph,
        params.kourami_ref,
        params.kourami_database,
        params.trimmer,
        params.adapter_fasta,
        params.save_trimmed_fail,
        params.save_merged
    )

    } else {
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
}
