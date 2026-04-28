process BWA_MEM_ALIGN_ALT_POSTALT {
    tag "$meta.sample"
    publishDir "$params.outdir/bwa-aln-postalt"
    label "bwa_mem_container"

    input:
    path reference
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_postalt.bam"), emit: bamfile_postalt
    path("*_postalt.bam.flagstat")
    path("*_postalt.bam.md5")

    script:
    """
    bwa mem -K 10000000 -t ${task.cpus} hs38DH.fa ${reads} > "${meta.sample}.bwamem.sam"
    k8 /usr/local/bin/bwa-postalt.js \
    hs38DH.fa.alt \
    ${meta.sample}.bwamem.sam > ${meta.sample}_postalt.sam
    rm "${meta.sample}.bwamem.sam"
    samtools view -bh -o ${meta.sample}_postalt.bam ${meta.sample}_postalt.sam
    rm "${meta.sample}_postalt.sam"
    samtools flagstat ${meta.sample}_postalt.bam > ${meta.sample}_postalt.bam.flagstat
    md5sum ${meta.sample}_postalt.bam > ${meta.sample}_postalt.bam.md5
    """
}
