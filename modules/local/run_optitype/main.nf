process RUN_OPTITYPE{
    tag "$meta.sample"

    publishDir "${params.outdir}/optitype_calls/${meta.sample}", mode: 'copy'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/optitype:1.3.5--0' :
        'quay.io/biocontainers/optitype:1.3.5--0' }"
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("optitype_calls"), emit: optitype_call // will be a directory called optitype_calls containing the tsv and pdf
    script:
    """
    gunzip *.gz
    OptiTypePipeline.py --input *.1.fq *.2.fq --verbose --dna --outdir ${meta.sample}
    mkdir -p optitype_calls
    cp "${meta.sample}"/*/*.tsv optitype_calls
    cp "${meta.sample}"/*/*.pdf optitype_calls
    rm *.fq
    """
}
