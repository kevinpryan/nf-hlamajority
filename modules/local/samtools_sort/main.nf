process samtools_sort{
    tag "$meta.sample"

    label "samtools_container"

    publishDir "$params.outdir/sort"

    input:
    tuple val(meta), path(bamfile)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: sortedbam

    script:
    """
    samtools sort -o ${meta.sample}.sorted.bam ${bamfile}
    """
}

/*
process samtools_sort{
    tag "$meta.sample"

    label 'samtools_container'
    publishDir "$params.outdir/sort"

    input:
    tuple val(meta), path(aln)
    path fasta

    output:
    tuple val(meta), path("*.sorted.{bam,cram}"), emit: sortedAln

    script:
    def cram_bam = (aln.extension == 'cram') ? "cram" : "bam"
    def ref_flag = (aln.extension == 'cram') ? "--reference ${fasta}" : ""

    """
    samtools sort \
             -o ${meta.sample}.sorted.${cram_bam} \
             -@ ${task.cpus} \
              ${ref_flag} \
             ${aln}
    """
}
*/

