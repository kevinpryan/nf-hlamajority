/*
process bam2fastq{
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
*/
/*
sambamba view \
    -f "bam" -h -p -l 0 -t ${task.cpus} \
    *.bam | sambamba sort -p -n -t ${task.cpus} -o - /dev/stdin |
    samtools fastq /dev/stdin \
    -1 "${meta.sample}_subset.1.fq" \
    -2 "${meta.sample}_subset.2.fq" \
    -0 /dev/null -s /dev/null -n
*/

process bam2fastq {
    tag "$meta.sample"
    label 'samtools_container'
    publishDir "$params.outdir/bam2fastq"

    input:
    tuple val(meta), path(reads), path(index) // reads is the bam file

    output:
    // This pattern matches both "sample.1.fq.gz" (PE) and "sample.fq.gz" (SE)
    tuple val(meta), path("*.fq.gz"), emit: convertedfastqs

    script:
    """
    mkdir -p tmp
    
    # Check BAM header/flags to see if it contains paired reads (Flag 1 = paired)
    # We check the first 1000 reads to be safe and fast.
    if samtools view -f 1 $reads | head -n 1 | grep -q .; then
        echo "Detected Paired-End BAM"
        MODE="PE"
    else
        echo "Detected Single-End BAM"
        MODE="SE"
    fi

    if [ "\$MODE" == "PE" ]; then
        samtools sort -l 9 -T tmp -@ ${task.cpus} -n ${reads} -o - | \\
        samtools fastq -c 6 /dev/stdin \\
            -1 "${meta.sample}.1.fq.gz" \\
            -2 "${meta.sample}.2.fq.gz" \\
            -0 /dev/null -s /dev/null -n
    else
        samtools sort -l 9 -T tmp -@ ${task.cpus} -n ${reads} -o - | \\
        samtools fastq -c 6 /dev/stdin \\
            -0 "${meta.sample}.fq.gz" -n
    fi
    """
}
