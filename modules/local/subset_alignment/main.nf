process SUBSET_ALIGNMENT {
    tag "$meta.sample"
    label "samtools_container"

    input:
    tuple val(meta), path(aln), path(index)
    path fasta

    output:
    tuple val(meta), path("${meta.sample}_subset.bam"), emit: subset_bam

    script:
    def ref_flag = (aln.extension == 'cram') ? "--reference ${fasta}" : ""
    """
    REGIONS=\$(samtools idxstats ${aln} | cut -f1 | grep -E '^(chr)?6(_.*)?\$|^HLA' | tr '\\n' ' ')
    samtools view -b \
             $ref_flag \
             -o ${meta.sample}_subset.bam \
             ${aln} \
             -@ ${task.cpus} \
             \$REGIONS '*'
    """
}

