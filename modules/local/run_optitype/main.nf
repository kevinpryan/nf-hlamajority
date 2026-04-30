process RUN_OPTITYPE {
    tag "$meta.sample"
    
    publishDir "${params.outdir}/optitype_calls/${meta.sample}", mode: 'copy'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/optitype:1.3.5--0' :
        'quay.io/biocontainers/optitype:1.3.5--0' }"

    input:
    tuple val(meta), path(reads) // Can be [read1, read2] or [read1]

    output:
    tuple val(meta), path("optitype_calls"), emit: optitype_call 

    script:
    """
    # Unzip the FASTQ files (OptiType works best with unzipped files)
    # This works for both SE (sample.fq.gz) and PE (sample.1.fq.gz / sample.2.fq.gz)
    gunzip -f *.gz
    OptiTypePipeline.py \
        --input *.fq \
        --dna \
        --verbose \
        --outdir ${meta.sample} \
        --prefix ${meta.sample}

    mkdir -p optitype_calls
    cp -r "${meta.sample}"/*.tsv optitype_calls/
    cp -r "${meta.sample}"/*.pdf optitype_calls/    
    rm *.fq
    """
}
