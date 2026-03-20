process RUN_KOURAMI_JAR{
    tag "$meta.sample"

    publishDir "${params.outdir}/kourami/${meta.sample}", mode: 'copy'
    input:
    tuple val(meta), path(bam_bai)
    path kourami_panel
    output:
    tuple val(meta), path("kourami_calls"), emit: kourami_result
    script:
    """
    mkdir -p kourami_calls
    java -jar /opt/wtsi-cgp/java/Kourami.jar --outfilePrefix ${meta.sample} -d ${kourami_panel} *.bam
    cp *.result kourami_calls
    """
}

process RUN_KOURAMI_PLACEHOLDER {
    tag "$meta.sample"
    publishDir "${params.outdir}/kourami/${meta.sample}", mode: 'copy'

    input:
    val meta

    output:
    tuple val(meta), path("kourami_calls"), emit: kourami_result

    script:
    """
    mkdir -p kourami_calls
    touch kourami_calls/${meta.sample}.result
    """
}
