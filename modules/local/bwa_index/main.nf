process bwa_index{
    tag "$meta.sample"
    publishDir "$params.outdir/bwa-aln"
    input:
    path reference
    //tuple val(meta), path(reads)
    //val reference_basename
    output:
    tuple val(meta), path("*.bam"), emit: bamfile
    path("*.bam.flagstat")
    script:
    """
    bwa index ${reference}
    """
}
