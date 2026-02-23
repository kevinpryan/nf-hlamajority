process BWA_INDEX {
    label 'bwa_mem_container'
    publishDir "${params.references_basedir}/${subdir}", mode: 'copy'

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
