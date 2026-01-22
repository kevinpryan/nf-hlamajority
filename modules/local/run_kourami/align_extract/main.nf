
process RUN_KOURAMI_ALIGN_EXTRACT{
    tag "$meta.sample"

    publishDir "${params.outdir}/kourami/${meta.sample}"
    input:
    tuple val(meta), path(bam), path(index)
    path database
    path reference
    output:
    tuple val(meta), path("*on_KouramiPanel.bam"), emit: kourami_alignment
    script:
    """
    bash alignAndExtract_hs38DH.sh -r hs38NoAltDH.fa -d "${database}" "${meta.sample}" ${bam}
    """
}
