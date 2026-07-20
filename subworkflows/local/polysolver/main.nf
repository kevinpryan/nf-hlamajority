//
// takes subset bam file from alt_align workflow and realigns without alt contigs as Polysolver does not accept alt alignments
// does not remove "chr" from header
// 

include { BAM_TO_FASTQ } from '../../../modules/local/bam_to_fastq'
include { BWA_REALIGN } from '../../../modules/local/bwa_realign'
include { SAMTOOLS_SORT_INDEX } from '../../../modules/local/samtools_sort_index'
include { RUN_POLYSOLVER } from '../../../modules/local/run_polysolver'
include { RUN_POLYSOLVER_PLACEHOLDER_SINGLE_END } from '../../../modules/local/run_polysolver'

workflow POLYSOLVER{
    /*
    need to realign without alt contigs
    */
    take: 
    subsetbam
    reference
    fasta_cram
    ch_novoalign
    ch_novolicense

    main:
    BAM_TO_FASTQ(
        subsetbam
    )

    BWA_REALIGN(
        BAM_TO_FASTQ.out.convertedfastqs,
        reference
        ) 

    SAMTOOLS_SORT_INDEX(
        BWA_REALIGN.out.realignbam,
        fasta_cram
        )

    RUN_POLYSOLVER(
        SAMTOOLS_SORT_INDEX.out.sortedAln,
        ch_novoalign,
        ch_novolicense
    )

    // Expected samples: keep only key
    SAMTOOLS_SORT_INDEX.out.sortedAln
        .map { meta, bam, index ->
            [meta.sample, meta]
        }
        .set { expected_samples }


    // Successful Polysolver outputs
    RUN_POLYSOLVER.out.polysolver_call
        .map { meta, result ->
            [meta.sample, result]
        }
        .set { successful_polysolver }


    // Find missing samples
    expected_samples
        .join(successful_polysolver, remainder:true)
        .filter { sample, meta, result ->
            result == null
        }
        .map { sample, meta, result ->
            meta
        }
        .branch {
            paired: it.single_end == false
            single: it.single_end == true
        }
        .set { failed_polysolver }
     
     RUN_POLYSOLVER_PLACEHOLDER_SINGLE_END(
        failed_polysolver.single
     )

    emit:
    calls = RUN_POLYSOLVER.out.polysolver_call
                              .mix(RUN_POLYSOLVER_PLACEHOLDER_SINGLE_END.out.polysolver_call)
    status = RUN_POLYSOLVER.out.run_status 
                               .mix(RUN_POLYSOLVER_PLACEHOLDER_SINGLE_END.out.run_status)
}
