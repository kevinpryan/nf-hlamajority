process samtools_sort_index{
    tag "$meta.sample"

    label 'samtools_container'
    publishDir "$params.outdir/sort_index"

    input:
    tuple val(meta), path(aln)
    path fasta

    output:
    tuple val(meta), path("*.sorted.{bam,cram}"), path("*.sorted.{bam.bai,cram.crai,bai,crai}"), emit: sortedAln

    script:
    def cram_bam = (aln.extension == 'cram') ? "cram" : "bam"
    def ref_flag = (aln.extension == 'cram') ? "--reference ${fasta}" : ""

    """
    samtools sort \
             -o ${meta.sample}.sorted.${cram_bam} \
             -@ ${task.cpus} \
              ${ref_flag} \
             ${aln}

    samtools index ${meta.sample}.sorted.${cram_bam}
    """
}
