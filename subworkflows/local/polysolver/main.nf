//
// takes subset bam file from alt_align workflow and realigns without alt contigs as Polysolver does not accept alt alignments
// does not remove "chr" from header
// 

include { BAM_TO_FASTQ } from '../../../modules/local/bam_to_fastq'
include { BWA_REALIGN } from '../../../modules/local/bwa_realign'
include { SAMTOOLS_SORT_INDEX } from '../../../modules/local/samtools_sort_index'
include { RUN_POLYSOLVER } from '../../../modules/local/run_polysolver'

workflow POLYSOLVER{
    /*
    need to realign without alt contigs
    */
    take: 
    subsetbam
    reference
    fasta_cram

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
        SAMTOOLS_SORT_INDEX.out.sortedAln
    )

    emit:
    calls = RUN_POLYSOLVER.out.polysolver_call
    status = RUN_POLYSOLVER.out.run_status 
}
