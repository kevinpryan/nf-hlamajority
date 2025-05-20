process realignwithoutAlt{
    publishDir "$params.outdir/realignwithoutAlt"
    
    input:
    tuple val(meta), path(reads)
    path reference
    val reference_basename

    output:
    tuple val(meta), path("*.bam"), emit: realignbam
    path("*_realign.sam.header")
    path("*_realign.sam.flagstat")

    script:
    """
    bwa mem -t ${task.cpus} -j ${reference_basename}.fa ${reads} > ${meta.sample}_realign.sam
    samtools view -H ${meta.sample}_realign.sam > ${meta.sample}_realign.sam.header
    samtools flagstat ${meta.sample}_realign.sam > ${meta.sample}_realign.sam.flagstat
    samtools view -bh -o ${meta.sample}.bam ${meta.sample}_realign.sam
    rm "${meta.sample}_realign.sam"
    """
}

/*
output:
    tuple val(meta), path("*_realign.sam"), emit: realignbam
    path("*_realign.sam.header")
    path("*_realign.sam.flagstat")
*/
