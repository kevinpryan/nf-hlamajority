process realignwithoutAlt{
    label 'bwa_mem_container'
    publishDir "$params.outdir/realignwithoutAlt"
    
    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.bam"), emit: realignbam
    path("*_realign.sam.header")
    path("*_realign.sam.flagstat")

    script:
    """
    bwa mem -t ${task.cpus} -j hs38DH.fa ${reads} > ${meta.sample}_realign.sam
    samtools view -H ${meta.sample}_realign.sam > ${meta.sample}_realign.sam.header
    samtools flagstat ${meta.sample}_realign.sam > ${meta.sample}_realign.sam.flagstat
    samtools view -bh -o ${meta.sample}.bam ${meta.sample}_realign.sam
    rm "${meta.sample}_realign.sam"
    """
}
