/*
process bwa_mem_align_alt_postalt{
    publishDir "$params.outdir/bwa-aln-postalt"
    input:
    path reference
    tuple val(meta), path(reads)
    val reference_basename
    output:
    tuple val(meta), path("*_postalt.bam"), emit: bamfile_postalt
    path("*_postalt.bam.flagstat")
    script:
    """
    bwa mem -t ${task.cpus} ${reference_basename}.fa ${reads} > "${meta.sample}.bwamem.sam"
    k8-linux /usr/local/bin/bwa-0.7.15/bwa-postalt.js \
    ${reference_basename}.fa.alt \
    ${meta.sample}.bwamem.sam > ${meta.sample}_postalt.sam
    rm "${meta.sample}.bwamem.sam"
    samtools view -bh -o ${meta.sample}_postalt.bam ${meta.sample}_postalt.sam
    rm "${meta.sample}_postalt.sam"
    samtools flagstat ${meta.sample}_postalt.bam > ${meta.sample}_postalt.bam.flagstat
    """
}
*/
process bwa_mem_align_alt_postalt{
    publishDir "$params.outdir/bwa-aln-postalt"
    label "bwa_mem_container"
    input:
    path reference
    tuple val(meta), path(reads)
    val reference_basename
    output:
    tuple val(meta), path("*_postalt.bam"), emit: bamfile_postalt
    path("*_postalt.bam.flagstat")
    script:
    """
    bwa mem -t ${task.cpus} ${reference_basename}.fa ${reads} > "${meta.sample}.bwamem.sam"
    k8 /usr/local/bin/bwa-postalt.js \
    ${reference_basename}.fa.alt \
    ${meta.sample}.bwamem.sam > ${meta.sample}_postalt.sam
    rm "${meta.sample}.bwamem.sam"
    samtools view -bh -o ${meta.sample}_postalt.bam ${meta.sample}_postalt.sam
    rm "${meta.sample}_postalt.sam"
    samtools flagstat ${meta.sample}_postalt.bam > ${meta.sample}_postalt.bam.flagstat
    """
}

