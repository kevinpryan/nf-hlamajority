#!/usr/bin/env nextflow

include { HLATYPING } from "./workflows/hlatyping"
include { samtools_sort_index as samtools_sort_index_before_subset } from "./modules/local/samtools_sort_index"
include { samtools_sort_index as samtools_sort_index_after_subset } from "./modules/local/samtools_sort_index"
include { bam2fastq } from "./modules/local/bam2fastq"
include { SUBSET_ALIGNMENT } from "./modules/local/subset_alignment"

def cram_ref = params.cram_fasta ? file(params.cram_fasta, checkIfExists: true) : []
ch_fasta_cram = Channel.value(cram_ref)
reference_dir    = params.reference_dir    ?: "${params.references_basedir}/bwakit/hs38DH*"
hla_la_graph     = params.hla_la_graph     ?: "${params.references_basedir}/hla-la"
kourami_database = params.kourami_database ?: "${params.references_basedir}/kourami/db"
kourami_ref      = params.kourami_ref      ?: "${params.references_basedir}/kourami/resources/hs38NoAltDH.fa*"

workflow {
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
        samtools_sort_index_before_subset(
                            ch_alignment,
                            ch_fasta_cram
                            ) 
        
        SUBSET_ALIGNMENT(samtools_sort_index_before_subset.out.sortedAln, 
                         ch_fasta_cram
                         )

        samtools_sort_index_after_subset(
                            SUBSET_ALIGNMENT.out.subset_bam,
                            ch_fasta_cram
                            )

        bam2fastq(samtools_sort_index_after_subset.out.sortedAln)
        
        ch_fastq = bam2fastq.out.convertedfastqs
   
    } else {
    println "fastq input..."
    Channel.fromPath(params.samplesheet, checkIfExists: true)
    | splitCsv( header:true )
    | map { row ->
        meta = row.subMap('sample')
        [meta, [
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)]]
    }
    | set { ch_fastq }
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
        ch_fasta_cram
    )
}
