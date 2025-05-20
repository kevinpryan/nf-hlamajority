process markduplicates{
    publishDir "$params.outdir/markduplicates"

    input:
    tuple val(meta), path(sortedbam), path(sortedbam_index)

    output:
    tuple val(meta), path("*_sorted_mdup.bam*"), emit: markdupbam
    path("*.sorted.mdup.bam.flagstat")

    script:
    """
    bammarkduplicates I=${sortedbam} O=${meta.sample}_sorted_mdup.bam index=1 rmdup=0
    samtools flagstat ${meta.sample}_sorted_mdup.bam > ${meta.sample}.sorted.mdup.bam.flagstat
    """
}

// old version     tuple val(meta), path("*_sorted_mdup.bam"), path("*_sorted_mdup.bam.bai"), emit: markdupbam

