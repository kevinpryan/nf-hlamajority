process markduplicates{
    label "biobambam2_container"
    publishDir "$params.outdir/markduplicates"

    input:
    tuple val(meta), path(sortedbam), path(sortedbam_index)

    output:
    tuple val(meta), path("*_sorted_mdup.bam*"), emit: markdupbam

    script:
    """
    bammarkduplicates I=${sortedbam} O=${meta.sample}_sorted_mdup.bam index=1 rmdup=0
    """
}

// old version     tuple val(meta), path("*_sorted_mdup.bam"), path("*_sorted_mdup.bam.bai"), emit: markdupbam

