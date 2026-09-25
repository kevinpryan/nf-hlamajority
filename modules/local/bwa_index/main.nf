process BWA_INDEX {
    label 'bwa_mem_container'
    publishDir "${params.references_basedir}", mode: 'copy',
        saveAs: { file -> "${subdir}/${file}" }

    input:
    path reference
    val subdir

    output:
    path("*.{amb,ann,bwt,pac,sa}")

    script:
    """
    bwa index ${reference}
    """
}
