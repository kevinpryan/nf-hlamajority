process bam2fastq{
    publishDir "$params.outdir/bam2fastq"

    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("*.fq.gz"), emit: convertedfastqs

    script:
    """
    mkdir -p tmp
    samtools sort -l 9 -T tmp -@ ${task.cpus} -n *.bam -o - |
    samtools fastq -c 6 /dev/stdin \
    -1 "${meta.sample}.1.fq.gz" \
    -2 "${meta.sample}.2.fq.gz" \
    -0 /dev/null -s /dev/null -n
    """
}
/*
sambamba view \
    -f "bam" -h -p -l 0 -t ${task.cpus} \
    *.bam | sambamba sort -p -n -t ${task.cpus} -o - /dev/stdin |
    samtools fastq /dev/stdin \
    -1 "${meta.sample}_subset.1.fq" \
    -2 "${meta.sample}_subset.2.fq" \
    -0 /dev/null -s /dev/null -n
*/
