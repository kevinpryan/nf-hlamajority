process MOSDEPTH {

    tag "$meta.sample"
    publishDir "$params.outdir/mosdepth"
    label "mosdepth_container"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("${meta.sample}.regions.bed.gz"), emit: mosdepth_output

    script:
    """
    mosdepth --by ${bed} \
             --no-per-base \
             --mapq 0 \
             ${meta.sample} \
             ${bam}
    """
}


process MEAN_COVERAGE {

    tag "$meta.sample"
    publishDir "$params.outdir/mosdepth"
    label "r_basic_container"

    input:
    tuple val(meta), path(mosdepth_output)

    output:
    tuple val(meta), path("${meta.sample}-mean-depth-hla-classI-exons-2-3.tsv"), emit: mean_depth
    script:
    """
    calculate_mean_depth_mosdepth.R --input ${mosdepth_output} \
                                    --samplename ${meta.sample} \
                                    --output ${meta.sample}-mean-depth-hla-classI-exons-2-3.tsv
    """
}

