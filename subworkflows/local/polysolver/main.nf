//
// takes subset bam file from alt_align workflow and realigns without alt contigs as Polysolver does not accept alt alignments
// does not remove "chr" from header
// 

include { BAM_TO_FASTQ } from '../../../modules/local/bam_to_fastq'
include { BWA_REALIGN } from '../../../modules/local/bwa_realign'
include { SAMTOOLS_SORT_INDEX } from '../../../modules/local/samtools_sort_index'
include { RUN_POLYSOLVER } from '../../../modules/local/run_polysolver'
include { RUN_POLYSOLVER_PLACEHOLDER_SINGLE_END } from '../../../modules/local/run_polysolver'
include { RUN_POLYSOLVER_PLACEHOLDER_MISSING_NOVOALIGN } from '../../../modules/local/run_polysolver'

workflow POLYSOLVER{
    /*
    need to realign without alt contigs
    */
    take: 
    subsetbam
    reference
    fasta_cram

    main:

    // branch subsetbam here to give single and paired-end
    subsetbam
        .branch { meta, bam, index ->
               single: meta.single_end == true
               paired: meta.single_end == false
             }
        .set{ subsetbam_branched }

    // run placeholder process outputting Polysolver outputs (NA) and status file
    RUN_POLYSOLVER_PLACEHOLDER_SINGLE_END(
        subsetbam_branched.single
    )

    // proceed with Polysolver on paired-end samples

    BAM_TO_FASTQ(
        subsetbam_branched.paired
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
        SAMTOOLS_SORT_INDEX.out.sortedAln
    )

    emit:
    calls = RUN_POLYSOLVER.out.polysolver_call
                              .mix(RUN_POLYSOLVER_PLACEHOLDER_SINGLE_END.out.polysolver_call)
    status = RUN_POLYSOLVER.out.run_status 
                               .mix(RUN_POLYSOLVER_PLACEHOLDER_SINGLE_END.out.run_status)
}

workflow POLYSOLVER_PLACEHOLDER {
    /*
    run if Novoalign binary is not provided
    */

    take:
    subsetbam

    main:

    RUN_POLYSOLVER_PLACEHOLDER_MISSING_NOVOALIGN(
        subsetbam
    )

    emit:
    calls = RUN_POLYSOLVER_PLACEHOLDER_MISSING_NOVOALIGN.out.polysolver_call
    status = RUN_POLYSOLVER_PLACEHOLDER_MISSING_NOVOALIGN.out.run_status
}
