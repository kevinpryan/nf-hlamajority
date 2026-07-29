process RUN_KOURAMI_JAR{
    tag "$meta.sample"

    publishDir "${params.outdir}/kourami/${meta.sample}", mode: 'copy'
    input:
    tuple val(meta), path(bam_bai)
    path kourami_panel

    output:
    tuple val(meta), path("kourami_calls"), emit: kourami_result
    tuple val(meta), path("${meta.sample}.kourami.STATUS.txt"), emit: run_status

    script:
    """
    mkdir -p kourami_calls
    java -jar /opt/wtsi-cgp/java/Kourami.jar --outfilePrefix ${meta.sample} -d ${kourami_panel} *.bam
    cp *.result kourami_calls
    echo "${meta.sample}\tKourami\tSUCCESS" > "${meta.sample}.kourami.STATUS.txt"
    """
}

process RUN_KOURAMI_PLACEHOLDER {
    tag "$meta.sample"
    publishDir "${params.outdir}/kourami/${meta.sample}", mode: 'copy'

    input:
    val meta

    output:
    tuple val(meta), path("kourami_calls"), emit: kourami_result
    tuple val(meta), path("${meta.sample}.kourami.STATUS.txt"), emit: run_status

    script:
    """
    mkdir -p kourami_calls
    touch kourami_calls/${meta.sample}.result
    echo "${meta.sample}\tKourami\tTOOL_FAILURE" > "${meta.sample}.kourami.STATUS.txt"
    """
}
