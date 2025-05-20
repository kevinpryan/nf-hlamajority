process samtools_index{
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