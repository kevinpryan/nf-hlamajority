//
// takes subset bam file from alt_align workflow and realigns without alt contigs as Polysolver does not accept alt alignments
// also removes "chr" from bam header for Polysolver compatibility
// 

include { samtools_sort } from '../../../modules/local/samtools_sort'
include { bam2fastq } from '../../../modules/local/bam2fastq'
include { realignwithoutAlt } from '../../../modules/local/realignwithoutAlt'
include { reheaderChr } from '../../../modules/local/reheaderChr'
include { samtools_index } from '../../../modules/local/samtools_index'

workflow prepPolysolver{
    /*
    need to realign without alt contigs and remove "chr" from bam header
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
        bam2fastq.out.subsetfastq,
        reference,
        reference_basename
        ) 
    samtools_sort(
        realignwithoutAlt.out.realignbam
        )
    samtools_index(
        samtools_sort.out.sortedbam
    )
    samtools_sorted_index = samtools_sort.out.sortedbam.join(samtools_index.out.bam_indexed, by: 0)
    reheaderChr(
        samtools_sorted_index
    )
}
