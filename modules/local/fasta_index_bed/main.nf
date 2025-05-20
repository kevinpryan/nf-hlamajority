process fasta_index_bed{
    publishDir "$params.outdir/fasta_bed"

    input:
    path reference
    val chromosome
    val reference_basename
    output:
    path("*.bed"), emit: fasta_bed

    script:
    """
    samtools faidx ${reference_basename}.fa
    awk 'BEGIN {FS="\\t"}; {print \$1 FS "0" FS \$2}' ${reference_basename}.fa.fai > ${reference_basename}.full.fa.bed
    grep ${chromosome} ${reference_basename}.full.fa.bed > ${reference_basename}.fa.bed
    """
}
