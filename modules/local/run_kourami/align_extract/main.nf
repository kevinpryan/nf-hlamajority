
process RUN_KOURAMI_ALIGN_EXTRACT{
    publishDir "${params.outdir}/kourami/${meta.sample}"
    input:
    tuple val(meta), path(bam_bai)
    path database
    path reference
    output:
    tuple val(meta), path("*on_KouramiPanel.bam"), emit: kourami_alignment
    script:
    """
    bash alignAndExtract_hs38DH.sh -r hs38NoAltDH.fa -d "${database}" "${meta.sample}" *.bam
    """
}
// oldversion:     bash alignAndExtract_hs38DH.sh -r hs38NoAltDH.fa -d "${database}" "${meta.sample}" ${meta.sample}_subset.sorted.bam

// can specify kourami database (has to be correctly formatted) with -d
//     HLA-LA.pl --BAM *.bam --graph /usr/local/bin/HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT --workingDir . --sampleID ${meta.sample} --maxThreads ${task.cpus}

