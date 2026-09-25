process samtools_sort{
    tag "$meta.sample"

    label "samtools_container"

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
