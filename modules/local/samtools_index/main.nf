process samtools_index{
    tag "$meta.sample"

    label "samtools_container"
    publishDir "$params.outdir/index"

    input:
    tuple val(meta), path(sortedbam)

    output:
    tuple val(meta), path("*.sorted.bam.bai"), emit: bam_indexed

    script:
    """
    samtools index ${sortedbam}
    """
} 
/*
process samtools_index{
    tag "$meta.sample"

    label 'samtools_container'
    publishDir "$params.outdir/index"

    input:
    tuple val(meta), path(aln)
    path fasta

    output:
    tuple val(meta), path("*.sorted.{bam.bai,cram.crai,bai,crai}"), emit: idx

    script:
    def cram_bam = (aln.extension == 'cram') ? "cram" : "bam"
    """
    samtools index ${meta.sample}.sorted.${cram_bam}
    """
}
*/
