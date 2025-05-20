process SORT_RESULTS{
    publishDir "$params.outdir/combined_results", mode: 'copy'
    input:
    each path(outfile)
    output:
    path("*sorted.tsv")

    script:
    """ 
    awk 'NR==1{print; next} {print | "sort"}' $outfile >  ${outfile.baseName}_sorted.tsv
    """
}

