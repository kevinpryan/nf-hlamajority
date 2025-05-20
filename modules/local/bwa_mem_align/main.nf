process bwa_mem_align{
    publishDir "$params.outdir/bwa-aln"
    input:
    path reference
    tuple val(meta), path(reads)
    val reference_basename
    output:
    tuple val(meta), path("*.bam"), emit: bamfile
    path("*.bam.flagstat")
    script:
    """
    bwa mem -t ${task.cpus} -j ${reference_basename}.fa ${reads} > "${meta.sample}.bwamem.sam"
    samtools view -bh -o ${meta.sample}.bam ${meta.sample}.bwamem.sam
    rm "${meta.sample}.bwamem.sam"
    samtools flagstat ${meta.sample}.bam > ${meta.sample}.bam.flagstat
    """
}
