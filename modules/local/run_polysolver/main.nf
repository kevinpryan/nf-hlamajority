process RUN_POLYSOLVER {
    tag "$meta.sample"
    publishDir "${params.outdir}/polysolver_calls/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(idx)
    path novoalign
    path novolicense

    output:
    tuple val(meta), path("polysolver_calls"), emit: polysolver_call
    tuple val(meta), path("counts*"), optional: true
    path("check.status.out.txt"), optional: true
    //path("${meta.sample}.polysolver_status.tsv"), emit: run_status
    tuple val(meta), path("${meta.sample}.polysolver.STATUS.txt"), emit: run_status

    when:
    !meta.single_end

    script:
    """
    export NOVOALIGN_DIR=\$PWD
    chmod +x novoalign

    if [[ -n "${novolicense}" && -s "${novolicense}" ]]; then
        echo "Using Novoalign license ${novolicense}"
    else
        echo "No Novoalign license provided"
    fi

    mkdir -p polysolver_calls
    mkdir -p tempdir

    if /home/polysolver/scripts/shell_call_hla_type ${bam} Unknown 0 hg38 FASTQ 0 ./tempdir; then
        mv tempdir/winners.hla.nofreq.txt polysolver_calls/
        printf "${meta.sample}\\tPASS\\n" > ${meta.sample}.polysolver_status.tsv

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
    val meta

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
