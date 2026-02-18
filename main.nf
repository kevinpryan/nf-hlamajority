#!/usr/bin/env nextflow

include { HLATYPING } from "./workflows/hlatyping"
include { SAMTOOLS_SORT_INDEX as SAMTOOLS_SORT_INDEX_BEFORE_INDEX } from "./modules/local/samtools_sort_index"
include { SAMTOOLS_SORT_INDEX as SAMTOOLS_SORT_INDEX_AFTER_INDEX } from "./modules/local/samtools_sort_index"
include { BAM_TO_FASTQ } from "./modules/local/bam_to_fastq"
include { SUBSET_ALIGNMENT } from "./modules/local/subset_alignment"


workflow {
    ch_fasta_cram = params.cram_fasta ? Channel.value(file(params.cram_fasta)) : Channel.value([])
    references_basedir = params.references_basedir ?: "${projectDir}/references"
    reference_dir = "${references_basedir}/bwakit/hs38DH*"
    hla_la_graph = "${references_basedir}/hla-la"
    kourami_database = "${references_basedir}/kourami/db"
    kourami_ref = "${references_basedir}/kourami/resources/hs38NoAltDH.fa*"
    weights = params.weights ?: "${projectDir}/assets/benchmarking_results_claeys_cleaned.csv"

    if (params.aligned) {
        println "params.aligned specified..."
        // --- ALIGNMENT BRANCH (BAM/CRAM) ---
        Channel.fromPath(params.samplesheet, checkIfExists: true)
        | splitCsv(header: true)
        | map { row ->
            def meta = row.subMap(['sample'])
            def alignment_file = file(row.aln, checkIfExists: true)
            
            // Basic validation: if CRAM is used, ensure --fasta was provided
            if (alignment_file.extension == 'cram' && !params.cram_fasta) {
                error "ERROR: CRAM file detected [${alignment_file.name}], but no reference FASTA provided via --fasta"
            }
            return [ meta, alignment_file ]
        }
        | set { ch_alignment }
        SAMTOOLS_SORT_INDEX_BEFORE_INDEX(
                            ch_alignment,
                            ch_fasta_cram
                            ) 
        
        SUBSET_ALIGNMENT(SAMTOOLS_SORT_INDEX_BEFORE_INDEX.out.sortedAln, 
                            ch_fasta_cram
                        )

        SAMTOOLS_SORT_INDEX_AFTER_INDEX(
                            SUBSET_ALIGNMENT.out.subset_bam,
                            ch_fasta_cram
                            )

        BAM_TO_FASTQ(SAMTOOLS_SORT_INDEX_AFTER_INDEX.out.sortedAln)
        
        ch_fastq = BAM_TO_FASTQ.out.convertedfastqs.map { meta, reads ->
            meta.single_end = !(reads instanceof List) || reads.size() == 1
            tuple(meta, reads)
        }
        ch_fastq.view()
    } else {
    println "fastq input..."
    Channel.fromPath(params.samplesheet, checkIfExists: true)
    | splitCsv( header:true, strip:true )
    | flatMap { row ->
        if (!row.sample || !row.fastq_1) {
                // Return empty list to skip this row (ignores trailing empty lines)
                return [] 
        }
    def fastq_1 = file(row.fastq_1, checkIfExists: true)
    def reads = [ fastq_1 ]

    if (row.fastq_2) {
        reads << file(row.fastq_2, checkIfExists: true)
    }

    def meta = row.subMap('sample')
    meta.single_end = (reads.size() == 1)

    return [ [ meta, reads ] ]
    }
    | set { ch_fastq }
    ch_fastq.view()
    }
// example ch_fastq: [[sample:3532, seq_type:dna], [/data4/kryan/misc/useful/nextflow/nf-hlatyping/testdir/gen_testdata/3532_subset_10000.1.fq.gz, /data4/kryan/misc/useful/nextflow/nf-hlatyping/testdir/gen_testdata/3532_subset_10000.2.fq.gz]]

    HLATYPING(
        ch_fastq,
        reference_dir,
        hla_la_graph,
        kourami_ref,
        kourami_database,
        params.trimmer,
        params.adapter_fasta,
        params.save_trimmed_fail,
        params.save_merged,
        ch_fasta_cram,
        weights,
        params.voting_method
    )
}
