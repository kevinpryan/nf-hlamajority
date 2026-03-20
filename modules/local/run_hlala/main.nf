process RUN_HLALA {
    tag "$meta.sample"

    publishDir "${params.outdir}/hlala_calls/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(index)
    path graphdir

    output:
    tuple val(meta), path("hlala_calls"), emit: hlala_call

    when:
    !meta.single_end

    script:
    """
    mkdir -p hlala_calls

    HLA-LA.pl \
        --BAM ${bam} \
        --customGraphDir ${graphdir} \
        --graph PRG_MHC_GRCh38_withIMGT \
        --workingDir . \
        --sampleID ${meta.sample} \
        --maxThreads ${task.cpus}

    cp "${meta.sample}/hla/R1_bestguess_G.txt" hlala_calls/
    """
}

process RUN_HLALA_PLACEHOLDER {
    tag "$meta.sample"

    publishDir "${params.outdir}/hlala_calls/${meta.sample}", mode: 'copy'

    input:
    val meta

    output:
    tuple val(meta), path("hlala_calls"), emit: hlala_call

    script:
    """
    mkdir -p hlala_calls
    echo -e "Locus\\tChromosome\\tAllele\\tQ1\\tQ2\\tAverageCoverage\\tCoverageFirstDecile\\tMinimumCoverage\\tproportionkMersCovered\\tLocusAvgColumnError\\tNColumns_UnaccountedAllele_fGT0.2\\tperfectG" > hlala_calls/R1_bestguess_G.txt
        
    for locus in A A B B C C DQA1 DQA1 DQB1 DQB1 DRB1 DRB1 DPA1 DPA1 DPB1 DPB1 DRB3 DRB3 DRB4 DRB4 E E F F G G; do
        echo -e "\${locus}\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA" >> hlala_calls/R1_bestguess_G.txt
    done
    """
}
