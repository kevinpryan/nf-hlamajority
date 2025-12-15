process bwa_mem_align_alt_postalt{
    tag "$meta.sample"
    publishDir "$params.outdir/bwa-aln-postalt"
    label "bwa_mem_container"
    input:
    path reference
    tuple val(meta), path(reads)
    //val reference_basename
    output:
    tuple val(meta), path("*_postalt.bam"), emit: bamfile_postalt
    path("*_postalt.bam.flagstat")
    script:
    """
    bwa mem -t ${task.cpus} hs38DH.fa ${reads} > "${meta.sample}.bwamem.sam"
    k8 /usr/local/bin/bwa-postalt.js \
    hs38DH.fa.alt \
    ${meta.sample}.bwamem.sam > ${meta.sample}_postalt.sam
    rm "${meta.sample}.bwamem.sam"
    samtools view -bh -o ${meta.sample}_postalt.bam ${meta.sample}_postalt.sam
    rm "${meta.sample}_postalt.sam"
    samtools flagstat ${meta.sample}_postalt.bam > ${meta.sample}_postalt.bam.flagstat
    """
}

