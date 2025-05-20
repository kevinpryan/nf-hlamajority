process samtools_sort_index{
    publishDir "$params.outdir/sort_index"

    input:
    tuple val(meta), path(bamfile)

    output:
    tuple val(meta), path("*.sorted.bam*"), emit: sortedbam

    script:
    """
    samtools sort -o ${meta.sample}.sorted.bam ${bamfile}
    samtools index ${meta.sample}.sorted.bam
    """
}
