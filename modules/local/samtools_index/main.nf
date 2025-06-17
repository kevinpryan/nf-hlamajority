process samtools_index{
    label "samtools_container"
    publishDir "$params.outdir/index"

    input:
    tuple val(meta), path(sortedbam)

    output:
    tuple val(meta), path("*.sorted.bam.bai"), emit: bam_indexed

    script:
    """
    samtools index ${sortedbam}
    """
}   
