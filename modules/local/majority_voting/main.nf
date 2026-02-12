process MAJORITY_VOTE{

    tag "$meta.sample"
    publishDir "$params.outdir/results"
    label "r_basic_container"

    input:
    tuple val(meta), path(outputs)
    path weights
    path depth
    val method

    output:
    tuple val(meta), path("*_all_calls_mhci.tsv"), emit: all_calls
    tuple val(meta), path("*_votes_mhci.tsv"), emit: votes
    tuple val(meta), path("*_votes_mhci_stats.tsv"), emit: vote_stats

    script:
    """
    parse_outputs_majority_vote.R --samplename ${meta.sample} \
                                  --optitype optitype_calls \
                                  --polysolver polysolver_calls \
                                  --hlala hlala_calls \
                                  --kourami kourami_calls \
                                  --weights ${weights} \
                                  --method $method \
                                  --depth $depth
    """
}
