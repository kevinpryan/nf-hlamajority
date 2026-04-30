#!/usr/bin/env nextflow

include { ALT_ALIGN } from "../../subworkflows/local/alt_align"
include { OPTITYPE } from "../../subworkflows/local/optitype"
include { POLYSOLVER } from "../../subworkflows/local/polysolver"
include { HLA_LA } from "../../subworkflows/local/hlala"
include { KOURAMI } from "../../subworkflows/local/kourami"
include { FASTP } from "../../modules/nf-core/fastp"
include { MAJORITY_VOTE } from "../../modules/local/majority_voting"
include { SORT_RESULTS } from "../../modules/local/sort_results"
include { MOSDEPTH } from "../../modules/local/mosdepth"
include { MEAN_COVERAGE } from "../../modules/local/mosdepth"

workflow HLATYPING {
    take:
    ch_fastq
    reference_dir
    hla_la_graph
    kourami_ref
    kourami_database
    trim
    adapter_fasta
    save_trimmed_fail
    save_merged
    fasta_cram
    weights
    voting_method
    reference_polysolver
    main:

    ref = file(reference_dir, checkIfExists: true)
    graph = file(hla_la_graph, checkIfExists: true)
    ref_kourami = file(kourami_ref, checkIfExists: true)
    db_kourami = file(kourami_database, checkIfExists: true)
    weights = file(weights, checkIfExists: true)
    mosdepth_bed = file("$projectDir/assets/hla-a-b-c-exons-2-3.bed", checkIfExists: true)
    method = params.voting_method
    ref_polysolver = file(reference_polysolver, checkIfExists: true)
    if (trim == true) {
        FASTP (
        ch_fastq,
        [],
        save_trimmed_fail,
        save_merged
        )

        ch_fastq_align = FASTP.out.reads
    } else {
        println "skipping trimming..."
        ch_fastq_align = ch_fastq
    }

    ALT_ALIGN(
        ch_fastq_align,
        ref,
        fasta_cram
    )
    
    MOSDEPTH(
        ALT_ALIGN.out,
        mosdepth_bed
    )
    
    MEAN_COVERAGE(
        MOSDEPTH.out.mosdepth_output
    )

    OPTITYPE(
        ALT_ALIGN.out
    )
    
    POLYSOLVER(
        ALT_ALIGN.out,
        ref_polysolver,
        fasta_cram
    )
    
    HLA_LA(
        ALT_ALIGN.out,
        graph
    )

    KOURAMI(
        ALT_ALIGN.out,
        db_kourami,
        ref_kourami
    )

    OPTITYPE.out.mix(KOURAMI.out, POLYSOLVER.out.calls, HLA_LA.out)
                .groupTuple(by: 0, size: 4)
                .set{ ch_hlatyping_outputs }

    ch_hlatyping_outputs
                    .map{meta, results ->
                        [ meta, results.collect { it } ]
                    }
                    .set{ ch_hlatyping_outputs_grouped }

    ch_hlatyping_outputs_grouped
    .join(MEAN_COVERAGE.out.mean_depth) 
    .set { ch_majority_vote_inputs }

    MAJORITY_VOTE(
        ch_majority_vote_inputs,
        weights,
        voting_method
    )

    MAJORITY_VOTE.out.votes.collectFile(name: 'nf_hlamajority_votes_combined.tsv', keepHeader: true, skip: 1, sort: { it[0] }) { it[1] }
                            .set{ votes_ch }
    
    MAJORITY_VOTE.out.vote_stats.collectFile(name: 'nf_hlamajority_stats_combined.tsv', keepHeader: true, skip: 1, sort: { it[0] }) { it[1] }
                                .set{ votes_stats_ch }
 
    MAJORITY_VOTE.out.all_calls.collectFile(name: 'nf_hlamajority_all_calls.tsv', keepHeader: true, skip: 1, sort: { it[0] }) { it[1] }
                                .set{ all_calls_ch }


    MEAN_COVERAGE.out.mean_depth.collectFile(name: 'nf_hlamajority_depth.tsv', keepHeader: true, skip: 1,  sort: { it[0] }) { it[1] }
                                .set{ mosdepth_combined_ch }

    mixed_ch = votes_ch.mix(votes_stats_ch, all_calls_ch, mosdepth_combined_ch)

    SORT_RESULTS(mixed_ch)
}
