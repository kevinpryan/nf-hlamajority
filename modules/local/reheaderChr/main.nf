process reheaderChr{
    publishDir "$params.outdir/reheaderChr", mode: 'copy'

    input:
    tuple val(meta), path(bamfile), path(bamfile_idx)

    output:
    //tuple val(meta), path("*_reheader.bam"), path("*_reheader.bam.bai"), path("*_reheader.bam.flagstat"), path("*_reheader.bam.header"), emit: reheaderbam
    tuple val(meta), path("*_reheader.bam"), path("*_reheader.bam.bai"), emit: reheaderbam
    path("*_reheader.bam.flagstat")
    path("*_reheader.bam.header")
    script:
    """
    bash reheader_chr.sh ${bamfile} ${meta.sample}_reheader.bam
    samtools flagstat ${meta.sample}_reheader.bam > ${meta.sample}_reheader.bam.flagstat
    samtools view -H ${meta.sample}_reheader.bam > ${meta.sample}_reheader.bam.header
    """
}
