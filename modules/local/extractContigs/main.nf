process extractContigs {
    publishDir "$params.outdir/extractContigs"

    input:
    path hlatypes
    path reference
    val reference_basename
    val chromosome
    output:
    path("bwakit-alt_contigs.txt"), emit: alt_contigs
    path("bwakit-hla_contigs.txt"), emit: hla_contigs
    script:
    """
    grep ${chromosome} ${reference_basename}.fa.alt | awk -F '\\t' '{print \$1}' > bwakit-alt_contigs.txt
    awk -F '\\t' '{print \$3}' ${hlatypes} > bwakit-hla_contigs.txt
    """
}
