process RUN_HLALA{
    publishDir "${params.outdir}/hlala_calls/${meta.sample}", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    path graphdir
    output:
    tuple val(meta), path("hlala_calls"), emit: hlala_call
    script:
    """
    mkdir -p hlala_calls
    HLA-LA.pl --BAM *sorted*.bam --customGraphDir ${graphdir} --graph PRG_MHC_GRCh38_withIMGT --workingDir . --sampleID ${meta.sample} --maxThreads ${task.cpus}
    cp "${meta.sample}/hla/R1_bestguess_G.txt" hlala_calls
    """
}
//     HLA-LA.pl --BAM *.bam --graph /usr/local/bin/HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT --workingDir . --sampleID ${meta.sample} --maxThreads ${task.cpus}

