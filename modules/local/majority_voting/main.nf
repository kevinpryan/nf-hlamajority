process MAJORITY_VOTE{
    publishDir "$params.outdir/results"

    input:
    tuple val(meta), path(outputs)
    path benchmark
    output:
    tuple val(meta), path("*_all_calls_mhci.tsv"), emit: all_calls
    tuple val(meta), path("*_majority_vote_mhci.tsv"), emit: majority_vote
    path("dummy_out8.txt")
    script:
    """
    Rscript ${params.rundir}/bin/parse_outputs_majority_vote.R --samplename ${meta.sample} --optitype optitype_calls --polysolver polysolver_calls --hlala hlala_calls --kourami kourami_calls --benchmark ${benchmark}
    """
}
