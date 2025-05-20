process samtools_sort{
    publishDir "$params.outdir/sort"

    input:
    tuple val(meta), path(bamfile)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: sortedbam

    script:
    """
    samtools sort -o ${meta.sample}.sorted.bam ${bamfile}
    """
}