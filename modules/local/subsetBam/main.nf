process subsetBam{
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
/*
process subsetBam2{
    publishDir "$params.outdir/subsetBam2"

    input:
    //tuple val(meta), path(bamfile), path(bamfileIndex)
    tuple val(meta), path(bam_bai)
    //path alt_chr6_contigs
    //path hla_contigs
    //path reference
    //path fasta_bed
    //path subset_regions
    output:
    tuple val(meta), path("*_subset.sorted.bam*"), emit: subsetbam
    path("*_subset.sorted.bam.flagstat")
    path("*_subset.sorted.bam.header")
    shell:
    '''

    echo 'subsetting regions...'
    # Get the list of chromosome 6 and its alt contigs
    chr6_contigs=$(samtools idxstats *.bam | awk '$1 == "chr6" || $1 ~ /^chr6_/ || $1 ~ /^HLA-/ {print $1}' | tr $'\n' ' ')
    #chr6_contigs=$(samtools idxstats *.bam | awk '$1 == "chr6" || $1 ~ /^chr6_/ || $1 ~ /^HLA-/ {print $1}' | tr '\n' ' ')    
    # Extract reads mapped to chromosome 6 and its alt contigs
    #samtools view -b -@ ~{num_threads} -o output_subset.bam ~{inputBAM} $chr6_contigs
    samtools view -b -@ -o output_subset.bam *.bam ${chr6_contigs}
    samtools sort -o !{meta.sample}_subset.sorted.bam output_subset.bam
    echo 'indexing sorted subset bam...'
    sambamba index !{meta.sample}_subset.sorted.bam
    echo 'viewing header...'
    samtools view -H !{meta.sample}_subset.sorted.bam > !{meta.sample}_subset.sorted.bam.header
    echo 'flagstat'
    samtools flagstat !{meta.sample}_subset.sorted.bam > !{meta.sample}_subset.sorted.bam.flagstat
    '''

}
*/

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
    chr6_contigs=\$(samtools idxstats "${bamfile}" | awk '/^chr6|^HLA-/{print \$1}' | tr \$'\\n' ' ')

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
/*
process subsetBam2 {
    tag "${meta.id}"
    label 'samtools_container'
    publishDir "$params.outdir/subsetBam2", mode: 'copy'

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
    // Ensure meta.sample exists or provide a fallback for prefix
    prefix = meta.sample ? "${meta.sample}_chr6HlaUnmapped_subset.sorted" : "${meta.id}_chr6HlaUnmapped_subset.sorted"

    """
    echo "Input BAM: ${bamfile}"
    echo "Output prefix: ${prefix}"
    echo "CPUs for task: ${task.cpus}"

    echo 'Identifying chr6, alt, and HLA contigs...'
    mapped_contigs=\$(samtools idxstats "${bamfile}" | awk '/^chr6|^HLA-/{print \$1}' | paste -s -d ' ')

    echo "Mapped contigs to subset (if any): '\${mapped_contigs}'"

    if [ -n "\${mapped_contigs}" ]; then
        echo "Extracting reads from: \${mapped_contigs} AND all unmapped reads."
        (samtools view -h -@ ${task.cpus - 1 < 1 ? 1 : task.cpus - 1} "${bamfile}" \${mapped_contigs}; samtools view -@ ${task.cpus - 1 < 1 ? 1 : task.cpus - 1} "${bamfile}" -f 4) | samtools sort -@ ${task.cpus} -o "${prefix}.bam" -
    else
        echo "No chr6/HLA contigs found or specified. Extracting only unmapped reads."
        samtools view -h -f 4 -@ ${task.cpus - 1 < 1 ? 1 : task.cpus - 1} "${bamfile}" | samtools sort -@ ${task.cpus} -o "${prefix}.bam" -
    fi
    echo 'Indexing sorted subset BAM...'
    samtools index -@ ${task.cpus} "${prefix}.bam"

    echo 'Viewing header...'
    samtools view -H "${prefix}.bam" > "${prefix}.bam.header"

    echo 'Calculating flagstat...'
    samtools flagstat -@ ${task.cpus} "${prefix}.bam" > "${prefix}.bam.flagstat"
    
    echo 'Calculating idxstats...'
    samtools idxstats -@ ${task.cpus} "${prefix}.bam" > "${prefix}.bam.idxstats"
    echo "Process finished."
    """
}
*/

/*
"""
    echo "Input BAM: ${bamfile}"
    echo "Output prefix: ${prefix}"
    echo "CPUs for task: ${task.cpus}"

    echo 'Identifying chr6, alt, and HLA contigs...'
    # Get the list of chromosome 6, its alt contigs, and HLA contigs
    # We do NOT include '*' here because unmapped reads are handled by -f 4
    # Using paste -s -d ' ' is often cleaner than tr for this
    mapped_contigs=\$(samtools idxstats "${bamfile}" | awk '/^chr6|^HLA-/{print \$1}' | paste -s -d ' ')

    echo "Mapped contigs to subset (if any): '\${mapped_contigs}'"

    # We will stream mapped reads (with header) and unmapped reads (without header)
    # together and pipe to sort.
    # The subshell (...) concatenates the standard output of the commands within it.

    if [ -n "\${mapped_contigs}" ]; then
        echo "Extracting reads from: \${mapped_contigs} AND all unmapped reads."
        (samtools view -h -@ ${task.cpus - 1 < 1 ? 1 : task.cpus - 1} "${bamfile}" \${mapped_contigs}; samtools view -@ ${task.cpus - 1 < 1 ? 1 : task.cpus - 1} "${bamfile}" -f 4) | samtools sort -@ ${task.cpus} -o "${prefix}.bam" -
    else
        echo "No chr6/HLA contigs found or specified. Extracting only unmapped reads."
        # If no mapped_contigs, just get unmapped reads (with header) and sort
        samtools view -h -f 4 -@ ${task.cpus - 1 < 1 ? 1 : task.cpus - 1} "${bamfile}" | samtools sort -@ ${task.cpus} -o "${prefix}.bam" -
    fi
    # Note on threads: samtools view and sort can use multiple threads.
    # The -@ for samtools view is for decompression/compression threads.
    # The -@ for samtools sort is for sorting threads.
    # Distributing task.cpus: give some to view if sort is also threaded.
    # If only one command is threaded, it can take all task.cpus.
    # Here, sort is the most demanding, so it gets full task.cpus.
    # The view commands get task.cpus-1 (min 1) to leave one for sort's main thread if it matters,
    # or just use fewer threads for view if piping. Simpler: give main threads to sort.
    # A simpler threading strategy for the pipe:
    # (samtools view -h --threads ${task.cpus} "${bamfile}" \${mapped_contigs}; samtools view --threads ${task.cpus} "${bamfile}" -f 4) | samtools sort -@ ${task.cpus} -o "${prefix}.bam" -
    # However, samtools view's --threads is for compression/decompression, not number of regions.
    # The default behavior with pipes and multiple samtools instances sharing task.cpus is generally fine.
    # Let's use a common approach for thread allocation in pipes:
    # For the IF branch:
    # (samtools view -h -@ $((task.cpus > 1 ? task.cpus / 2 : 1)) "${bamfile}" \${mapped_contigs}; samtools view -@ $((task.cpus > 1 ? task.cpus / 2 : 1)) "${bamfile}" -f 4) | samtools sort -@ ${task.cpus} -m 2G -o "${prefix}.bam" -
    # For the ELSE branch:
    # samtools view -h -f 4 -@ $((task.cpus > 1 ? task.cpus / 2 : 1)) "${bamfile}" | samtools sort -@ ${task.cpus} -m 2G -o "${prefix}.bam" -
    # Let's revert to simpler threading for now, samtools is pretty good at managing internal threads with -@
    # If the mapped_contigs list is empty, the first samtools view in the 'if' block might error or output all reads.
    # The if/else structure handles this: if no mapped_contigs, only extract unmapped.

    echo 'Indexing sorted subset BAM...'
    samtools index -@ ${task.cpus} "${prefix}.bam"

    echo 'Viewing header...'
    samtools view -H "${prefix}.bam" > "${prefix}.bam.header"

    echo 'Calculating flagstat...'
    samtools flagstat -@ ${task.cpus} "${prefix}.bam" > "${prefix}.bam.flagstat"

    echo "Process finished."
    """
*/

//     samtools view -o ${meta.sample}_subset.bam -b *.bam chr6 chr6_ chrUn chrUn_ HLA "*"

// samtools view -o ${meta.sample}_subset.bam -b *.bam chr6:28509970-33480727 chr6* HLA* chrUn* "*"

/*
    samtools view -o ${meta.sample}_subset.bam -b ${bamfile} chr6 chr6_* HLA* "*"


    samtools view -o ${meta.sample}_subset.bam -b ${bamfile} -F 256 chr6:29500000-33200000 chr6_* HLA*

samtools view -o ${meta.sample}_subset1.bam -b ${bamfile} -L ${fasta_bed}
    samtools view -o ${meta.sample}_subset2.bam -b ${bamfile} "*"
    samtools view -o ${meta.sample}_subset3.bam -b ${bamfile} 'chr6:28509970-33480727'
    samtools merge ${meta.sample}_subset.bam ${meta.sample}_subset1.bam ${meta.sample}_subset2.bam ${meta.sample}_subset3.bam
    samtools sort -o ${meta.sample}_subset.sorted.bam ${meta.sample}_subset.bam
    sambamba index ${meta.sample}_subset.sorted.bam
    samtools view -H ${meta.sample}_subset.sorted.bam > ${meta.sample}_subset.sorted.bam.header
    samtools flagstat ${meta.sample}_subset.sorted.bam > ${meta.sample}_subset.sorted.bam.flagstat
*/
