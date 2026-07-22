process RUN_POLYSOLVER {
    tag "$meta.sample"
    publishDir "${params.outdir}/polysolver_calls/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(idx)

    output:
    tuple val(meta), path("polysolver_calls"), emit: polysolver_call
    tuple val(meta), path("counts*"), optional: true
    tuple val(meta), path("${meta.sample}.polysolver.STATUS.txt"), emit: run_status

    when:
    !meta.single_end

    script:
    """
    mkdir -p polysolver_calls
    mkdir -p tempdir

    if /home/polysolver/scripts/shell_call_hla_type ${bam} Unknown 0 hg38 FASTQ 0 ./tempdir; then
        mv tempdir/winners.hla.nofreq.txt polysolver_calls/
        printf "${meta.sample}\\tPolysolver\\tSUCCESS\\n" > "${meta.sample}.polysolver.STATUS.txt"

    else
        printf "HLA-A\\tNA\\tNA\\n" > polysolver_calls/winners.hla.nofreq.txt
        printf "HLA-B\\tNA\\tNA\\n" >> polysolver_calls/winners.hla.nofreq.txt
        printf "HLA-C\\tNA\\tNA\\n" >> polysolver_calls/winners.hla.nofreq.txt
        printf "${meta.sample}\\tPolysolver\\tTOOL_FAILURE\\n" > "${meta.sample}.polysolver.STATUS.txt"
    fi

    rm -rf tempdir
    """
}

process RUN_POLYSOLVER_PLACEHOLDER_SINGLE_END {
    tag "$meta.sample"

    publishDir "${params.outdir}/polysolver_calls/${meta.sample}", mode: 'copy'

    input:
    //val meta
    tuple val(meta), path(bam), path(idx)

    output:
    tuple val(meta), path("polysolver_calls"), emit: polysolver_call
    tuple val(meta), path("${meta.sample}.polysolver.STATUS.txt"), emit: run_status

    script:
    """
    mkdir -p polysolver_calls
    printf "HLA-A\\tNA\\tNA\\n" > polysolver_calls/winners.hla.nofreq.txt
    printf "HLA-B\\tNA\\tNA\\n" >> polysolver_calls/winners.hla.nofreq.txt
    printf "HLA-C\\tNA\\tNA\\n" >> polysolver_calls/winners.hla.nofreq.txt

    echo "${meta.sample}\tPolysolver\tSKIPPED_SINGLE_END" > "${meta.sample}.polysolver.STATUS.txt"
    """
}

process RUN_POLYSOLVER_PLACEHOLDER_MISSING_NOVOALIGN {
    tag "$meta.sample"

    publishDir "${params.outdir}/polysolver_calls/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(idx)

    output:
    tuple val(meta), path("polysolver_calls"), emit: polysolver_call
    tuple val(meta), path("${meta.sample}.polysolver.STATUS.txt"), emit: run_status

    script:
    """
    mkdir -p polysolver_calls
    printf "HLA-A\\tNA\\tNA\\n" > polysolver_calls/winners.hla.nofreq.txt
    printf "HLA-B\\tNA\\tNA\\n" >> polysolver_calls/winners.hla.nofreq.txt
    printf "HLA-C\\tNA\\tNA\\n" >> polysolver_calls/winners.hla.nofreq.txt

    echo "${meta.sample}\tPolysolver\tSKIPPED_MISSING_NOVOALIGN" > "${meta.sample}.polysolver.STATUS.txt"
    """
}

