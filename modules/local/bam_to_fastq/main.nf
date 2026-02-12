process BAM_TO_FASTQ {
    tag "$meta.sample"
    label 'samtools_container'
    publishDir "$params.outdir/bam2fastq"

    input:
    tuple val(meta), path(reads), path(index)
    output:
    tuple val(meta), path("*.fq.gz"), emit: convertedfastqs

    script:
    """
    mkdir -p tmp
    samtools sort -l 9 -T tmp -@ ${task.cpus} -n ${reads} -o - |
    samtools fastq -c 6 /dev/stdin \
    -1 "${meta.sample}.1.fq.gz" \
    -2 "${meta.sample}.2.fq.gz" \
    -0 /dev/null -s /dev/null -n
    """
}
