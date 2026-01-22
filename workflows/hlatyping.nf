#!/usr/bin/env nextflow

include {alt_align} from "../subworkflows/local/alt_align"
include { prepPolysolver } from "../subworkflows/local/prepPolysolver"
include { optitype } from "../subworkflows/local/optitype"
include { polysolver } from "../subworkflows/local/polysolver"
include { hlala } from "../subworkflows/local/hlala"
include { kourami } from "../subworkflows/local/kourami"
include { FASTP } from "../modules/nf-core/fastp"
include { MAJORITY_VOTE } from "../modules/local/majority_voting"
include { SORT_RESULTS } from "../modules/local/sort_results"
include { MOSDEPTH } from "../modules/local/mosdepth"
include { MEAN_COVERAGE } from "../modules/local/mosdepth"

workflow HLATYPING {
    take:
    ch_fastq
    reference_dir
    hla_la_graph
    kourami_ref
    kourami_database
    trimmer
    adapter_fasta
    save_trimmed_fail
    save_merged
    ch_fasta_cram

    main:
    ch_ref = file(reference_dir, checkIfExists: true)
    ch_graph = file(hla_la_graph, checkIfExists: true)
    ch_ref_kourami = file(kourami_ref, checkIfExists: true)
    ch_db_kourami = file(kourami_database, checkIfExists: true)
    ch_benchmark = file("$projectDir/assets/benchmarking_results_claeys.csv", checkIfExists: true)
    ch_mosdepth_bed = file("$projectDir/assets/hla-a-b-c-exons-2-3.bed", checkIfExists: true)

    if (trimmer == 'fastp') {
    FASTP (
    ch_fastq,
    [],
    save_trimmed_fail,
    save_merged
    )
    ch_fastq_align = FASTP.out.reads
    } else {
    ch_fastq_align = ch_fastq
    }
    alt_align(
        ch_fastq_align,
        ch_ref,
        ch_fasta_cram
    )

    mosdepth(
        alt_align.out,
        ch_mosdepth_bed
    )
    
    MEAN_COVERAGE(
        mosdepth.out.mosdepth_output
    )
    optitype(
        alt_align.out
    )
    polysolver(
        alt_align.out,
        ch_ref,
        ch_fasta_cram
    )
    hlala(
        alt_align.out,
        ch_graph
    )
     
    kourami(
        alt_align.out,
        ch_db_kourami,
        ch_ref_kourami
    )
    optitype.out.mix(kourami.out, polysolver.out, hlala.out)
           .groupTuple(by: 0, size: 4)
           .set{ ch_hlatyping_outputs }
    /*
    optitype.out.mix(kourami.out, polysolver.out)
           .groupTuple(by: 0, size: 3)
           .set{ ch_hlatyping_outputs }
    */
    ch_hlatyping_outputs
                    .map{meta, results ->
                        //[ meta, results.collect { it.getParent() } ]
                        [ meta, results.collect { it } ]
                    }
                    .set{ ch_hlatyping_outputs_grouped }
    ch_hlatyping_outputs_grouped.view()
    MAJORITY_VOTE(
        ch_hlatyping_outputs_grouped,
        ch_benchmark
    )
    MAJORITY_VOTE.out.majority_vote.collectFile(name: 'nf_hlamajority_results_majority_vote_combined.tsv', keepHeader: true, skip: 1, sort: { it[0] }) { it[1] }.set{ majority_ch }
    MAJORITY_VOTE.out.all_calls.collectFile(name: 'nf_hlamajority_hlatyping_results_all_calls.tsv', keepHeader: true, skip: 1, sort: { it[0] }) { it[1] }.set{all_calls_ch}
    mixed_ch = majority_ch.mix(all_calls_ch)
    SORT_RESULTS(mixed_ch)
    MEAN_COVERAGE.out.mean_depth.collectFile(name: 'nf_hlamajority_mean_depth_exons2_3_hla_classI.csv', keepHeader: true, skip: 1)
}
