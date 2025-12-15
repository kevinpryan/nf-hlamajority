process subsetBam{
    tag "$meta.sample"

    publishDir "$params.outdir/subsetBam"

    input:
    //tuple val(meta), path(bamfile), path(bamfileIndex)
    tuple val(meta), path(bam_bai) 
    path alt_chr6_contigs
    path hla_contigs
    path reference
    path fasta_bed
    path subset_regions
    output:
    tuple val(meta), path("*_subset.sorted.bam*"), emit: subsetbam
    path("*_subset.sorted.bam.flagstat")
    path("*_subset.sorted.bam.header") 
    shell:
    '''
    echo 'subsetting regions...'
    samtools view -h *.bam | \
    awk 'NR==FNR {a[$1]; next} /^@/ || $3 in a' !{subset_regions} - | \
    samtools view -Sb - > !{meta.sample}_subset1.bam 
    echo 'subsetting unmapped reads...'
    samtools view -o !{meta.sample}_subset2.bam -b !{meta.sample}_sorted_mdup.bam "*"
    echo 'merging bams...'
    samtools merge !{meta.sample}_subset.bam !{meta.sample}_subset1.bam !{meta.sample}_subset2.bam 
    echo 'sorting merged bams...'
    samtools sort -o !{meta.sample}_subset.sorted.bam !{meta.sample}_subset.bam
    echo 'indexing sorted merged bam...'
    sambamba index !{meta.sample}_subset.sorted.bam
    echo 'viewing header...'
    samtools view -H !{meta.sample}_subset.sorted.bam > !{meta.sample}_subset.sorted.bam.header
    echo 'flagstat'
    samtools flagstat !{meta.sample}_subset.sorted.bam > !{meta.sample}_subset.sorted.bam.flagstat
    '''
}

process subsetBam2{
    tag "${meta.sample}"
    publishDir "$params.outdir/subsetBam2", mode: 'copy'
    label 'samtools_container'

    cpus 4 // Example: request 4 cpus

    input:
    tuple val(meta), path(bam_bai) // bam_bai is a list: [bam_file, bai_file]

    output:
    tuple val(meta), path("${prefix}.bam*"), emit: subsetbam // Will catch .bam and .bam.bai
    path("${prefix}.bam.flagstat"), emit: flagstat
    path("${prefix}.bam.header"), emit: header
    path("${prefix}.bam.idxstats"), emit: idxstats

    script:
    def bamfile = bam_bai[0] // Get the BAM file from the input tuple
    prefix = "${meta.sample}_subset.sorted" // Use meta.sample for unique prefix

    """
    echo "Input BAM: ${bamfile}"
    echo "Output prefix: ${prefix}"
    echo "CPUs for task: ${task.cpus}"

    echo 'Subsetting regions...'
    chr6_contigs=\$(samtools idxstats "${bamfile}" | awk '/^chr6|^HLA-/' | cut -f 1 | tr \$'\\n' ' ')

    echo "Contigs to subset: '\${chr6_contigs}'"

    samtools view -b -@ ${task.cpus} -o output_temp_subset.bam "${bamfile}" \${chr6_contigs}

    echo 'Sorting subset BAM...'
    samtools sort -@ ${task.cpus} -o "${prefix}.bam" output_temp_subset.bam

    echo 'Indexing sorted subset BAM...'
    samtools index -@ ${task.cpus} "${prefix}.bam"

    echo 'Viewing header...'
    samtools view -H "${prefix}.bam" > "${prefix}.bam.header"

    echo 'Calculating flagstat...'
    samtools flagstat "${prefix}.bam" > "${prefix}.bam.flagstat"
    echo 'Calculating idxstats...'
    samtools idxstats -@ ${task.cpus} "${prefix}.bam" > "${prefix}.bam.idxstats"

    echo "Process finished."
    """
}
