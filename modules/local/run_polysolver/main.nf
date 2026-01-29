/*
process RUN_POLYSOLVER{
    tag "$meta.sample"

    publishDir "${params.outdir}/polysolver_calls/${meta.sample}", mode: 'copy'
    input:
    tuple val(meta), path(bam), path(idx)
    output:
    tuple val(meta), path("polysolver_calls"), emit: polysolver_call // file of interest is winners.hla.nofreq.txt
    tuple val(meta), path("counts*")
    path("check.status.out.txt")
    script:
    """
    mkdir -p polysolver_calls
    mkdir -p tempdir
    /home/polysolver/scripts/shell_call_hla_type ${bam} Unknown 0 hg38 ILMFQ 0 ./ tempdir
    mv winners.hla.nofreq.txt polysolver_calls
    rm -rf tempdir
    """
}
*/
process RUN_POLYSOLVER {
    tag "$meta.sample"
    publishDir "${params.outdir}/polysolver_calls/${meta.sample}", mode: 'copy', pattern: "winners*"

    input:
    tuple val(meta), path(bam), path(idx)

    output:
    tuple val(meta), path("polysolver_calls"), emit: polysolver_call
    tuple val(meta), path("counts*"), optional: true
    path("check.status.out.txt"), optional: true
    path("${meta.sample}.polysolver_status.tsv"), emit: run_status

    script:
    """
    mkdir -p polysolver_calls
    mkdir -p tempdir

    if /home/polysolver/scripts/shell_call_hla_type ${bam} Unknown 0 hg38 ILMFQ 0 ./ tempdir; then
        mv winners.hla.nofreq.txt polysolver_calls/
        printf "${meta.sample}\\tPASS\\n" > ${meta.sample}.polysolver_status.tsv

    else
        printf "HLA-A\\tNA\\tNA\\n" > polysolver_calls/winners.hla.nofreq.txt
        printf "HLA-B\\tNA\\tNA\\n" >> polysolver_calls/winners.hla.nofreq.txt
        printf "HLA-C\\tNA\\tNA\\n" >> polysolver_calls/winners.hla.nofreq.txt
        printf "${meta.sample}\\tFAIL\\n" > ${meta.sample}.polysolver_status.tsv
    fi

    rm -rf tempdir
    """
}
