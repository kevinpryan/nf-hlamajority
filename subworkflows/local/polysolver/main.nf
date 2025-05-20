//
// takes subset bam file from alt_align workflow and realigns without alt contigs as Polysolver does not accept alt alignments
// does not remove "chr" from header
// 

include { samtools_sort } from '../../../modules/local/samtools_sort'
include { bam2fastq } from '../../../modules/local/bam2fastq'
include { realignwithoutAlt } from '../../../modules/local/realignwithoutAlt'
include { samtools_index } from '../../../modules/local/samtools_index'
include { RUN_POLYSOLVER } from '../../../modules/local/run_polysolver'
include { samtools_sort_index } from '../../../modules/local/samtools_sort_index'
workflow polysolver{
    /*
    need to realign without alt contigs
    */
    take: 
    subsetbam
    reference
    reference_basename

    main:
    bam2fastq(
        subsetbam
    )
    realignwithoutAlt(
        bam2fastq.out.convertedfastqs,
        reference,
        reference_basename
        ) 
    samtools_sort_index(
        realignwithoutAlt.out.realignbam
        )
    RUN_POLYSOLVER(
        samtools_sort_index.out.sortedbam
    )
    emit:
    RUN_POLYSOLVER.out.polysolver_call 
}
