#!/usr/bin/env nextflow

params.imgt_version = "3.63.0"
params.imgt_commit  = "8382fbe"
params.build_references = false
params.kourami_commit = "545c770"
params.references_basedir = "references"
params.reference_dir = "${params.references_basedir}/bwakit/hs38DH*"
params.hla_la_graph = "${params.references_basedir}/hla-la"
params.kourami_database = "${params.references_basedir}/kourami/custom_db/3.63.0/"
params.kourami_ref = "${params.references_basedir}/kourami/resources/hs38NoAltDH.fa*"
params.trim = true
params.ref_polysolver = "${params.references_basedir}/polysolver/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna*"
params.novoalign = "${projectDir}/bin/novoalign"
params.novolicense = null

params.hs38noaltdh_fa_md5 = null
params.hs38dh_fa_md5      = null
params.polysolver_fna_md5 = null
params.hla_la_tar_md5 = null
//params.hla_la_prg_tar = 'PRG_MHC_GRCh38_withIMGT.tar.gz'
params.hla_la_prg_tar_md5 = null

// Eenforcement of novoalign placement
def expected_novoalign = file("${projectDir}/bin/novoalign").toString()
def provided_novoalign = file(params.novoalign).toString()

if (expected_novoalign != provided_novoalign) {
    exit 1, """
ERROR: Novoalign binary must be located exactly at:
    ${expected_novoalign}

You provided:
    ${provided_novoalign}
"""
}

if (!file(expected_novoalign).exists()) {
    exit 1, "ERROR: Novoalign binary not found at: ${expected_novoalign}"
}

// Enforcement of the optional novoalign license location
if (params.novolicense) {

    def expected_license = file("${projectDir}/bin/novoalign.lic").toString()
    def provided_license = file(params.novolicense).toString()

    if (expected_license != provided_license) {
        exit 1, """
ERROR: Novoalign license must be located exactly at:
    ${expected_license}

You provided:
    ${provided_license}
"""
    }

    if (!file(expected_license).exists()) {
        exit 1, "ERROR: Novoalign license not found at: ${expected_license}"
    }
}

log.info """
Novoalign license (optional for novoalign binaries v3 or less but required for v4+):
    If provided, it must be located at exactly:
        ${projectDir}/bin/novoalign.lic

Example usage:
    --novolicense ${projectDir}/bin/novoalign.lic

The Novoalign binary must be downloaded from www.novocraft.com by the user and should be located at:
        ${projectDir}/bin/novoalign

No other paths are accepted.
"""

include { REFERENCES } from "./workflows/references"
include { HLATYPING } from "./workflows/hlatyping"
include { SAMTOOLS_SORT_INDEX as SAMTOOLS_SORT_INDEX_BEFORE_INDEX } from "./modules/local/samtools_sort_index"
include { SAMTOOLS_SORT_INDEX as SAMTOOLS_SORT_INDEX_AFTER_INDEX } from "./modules/local/samtools_sort_index"
include { BAM_TO_FASTQ } from "./modules/local/bam_to_fastq"
include { SUBSET_ALIGNMENT } from "./modules/local/subset_alignment"

workflow {
    if (!params.outdir) {
        exit 1, "Pipeline parameter '--outdir' is mandatory. Please provide a path for the output directory."
    }
    references_basedir = params.references_basedir
    
    if (params.build_references){
        log.info "Mode: Building References (IMGT v${params.imgt_version})"
        
        def missing_md5 = []

         if (!params.hs38noaltdh_fa_md5) missing_md5 << '--hs38noaltdh_fa_md5'
         if (!params.hs38dh_fa_md5)      missing_md5 << '--hs38dh_fa_md5'
         if (!params.polysolver_fna_md5) missing_md5 << '--polysolver_fna_md5'
 
         // Only required when auto-downloading HLA-LA tarball
         if (!params.hla_la_prg_tar && !params.hla_la_tar_md5)
             missing_md5 << '--hla_la_tar_md5'
 
         if (missing_md5) {
             exit 1, "Missing required MD5 parameters: ${missing_md5.join(', ')}"
         }
 
         // Optional soft warning for custom tarball
         if (params.hla_la_prg_tar && !params.hla_la_prg_tar_md5) {
             log.warn "No --hla_la_prg_tar_md5 provided; skipping integrity check on custom tarball."
         }
        REFERENCES(
                  params.references_basedir,
                  params.imgt_commit,
                  params.imgt_version,
                  params.kourami_commit,
                  params.hs38noaltdh_fa_md5,
                  params.hs38dh_fa_md5,
                  params.hla_la_tar_md5,
                  params.hla_la_prg_tar_md5,
                  params.polysolver_fna_md5
              )
    } else {
        log.info "Mode: Running nf-hlamajority"
        ch_fasta_cram = params.cram_fasta ? Channel.value(file(params.cram_fasta)) : Channel.value([])
        
        reference_dir = params.reference_dir
        hla_la_graph = params.hla_la_graph
        kourami_database = params.kourami_database
        kourami_ref = params.kourami_ref
        weights = params.weights

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
        trim = false
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
    trim = params.trim

    // issue message if single-end data is provided
    ch_fastq
        .filter { meta, reads -> meta.single_end }
        .map { meta, reads -> meta.sample }
        .collect()
        .map { samples ->

            if (samples.size() > 0) {
                log.warn """
                Single-end data detected.

                ${samples.size()} sample(s) will only be processed by OptiType.
                Other HLA typing tools require paired-end reads.
                """
            }

            samples
        }

    }

   ch_novoalign    = Channel.value(file(params.novoalign))
   ch_novolicense = params.novolicense \
    ? Channel.value(file(params.novolicense)) \
    : Channel.value([])

// example ch_fastq: [[sample:3532, seq_type:dna], [/data4/kryan/misc/useful/nextflow/nf-hlatyping/testdir/gen_testdata/3532_subset_10000.1.fq.gz, /data4/kryan/misc/useful/nextflow/nf-hlatyping/testdir/gen_testdata/3532_subset_10000.2.fq.gz]]

    HLATYPING(
        ch_fastq,
        reference_dir,
        hla_la_graph,
        kourami_ref,
        kourami_database,
        trim,
        params.adapter_fasta,
        params.save_trimmed_fail,
        params.save_merged,
        ch_fasta_cram,
        weights,
        params.voting_method,
        params.ref_polysolver,
        ch_novoalign,
        ch_novolicense
    )
  }
}
