process RUN_POLYSOLVER{
    tag "$meta.sample"

    publishDir "${params.outdir}/polysolver_calls/${meta.sample}", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("polysolver_calls"), emit: polysolver_call // file of interest is winners.hla.nofreq.txt
    tuple val(meta), path("counts*")
    path("check.status.out.txt")
    script:
    """
    mkdir -p polysolver_calls
    mkdir -p tempdir
    /home/polysolver/scripts/shell_call_hla_type *.bam Unknown 0 hg38 ILMFQ 0 ./ tempdir
    mv winners.hla.nofreq.txt polysolver_calls
    rm -rf tempdir
    """
}
