process RUN_KOURAMI_ALIGN_EXTRACT{
    tag "$meta.sample"

    publishDir "${params.outdir}/kourami/${meta.sample}"

    input:
    tuple val(meta), path(bam), path(index)
    path database
    path reference

    output:
    tuple val(meta), path("*on_KouramiPanel.bam"), emit: kourami_alignment

    when:
    !meta.single_end

    script:
    """
    bash alignAndExtract_hs38DH.sh -r hs38NoAltDH.fa -d "${database}" "${meta.sample}" ${bam}
    """
}

process RUN_KOURAMI_PLACEHOLDER_SE {
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
    echo "${meta.sample}\tKourami\tSKIPPED_SINGLE_END" > "${meta.sample}.kourami.STATUS.txt"
    """
}
