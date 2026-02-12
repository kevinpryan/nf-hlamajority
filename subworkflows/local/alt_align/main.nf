//
// Alt-aware alignment with bwa-mem and postprocessing script
// Also subsets output bam file to relevant regions for HLA typing
//

include { BWA_MEM_ALIGN_ALT_POSTALT } from '../../../modules/local/bwa_mem_align_alt_postalt'
include { SAMTOOLS_SORT_INDEX as SAMTOOLS_SORT_INDEX_BEFORE_SUBSET } from '../../../modules/local/samtools_sort_index'
include { SAMTOOLS_SORT_INDEX as SAMTOOLS_SORT_INDEX_AFTER_SUBSET } from '../../../modules/local/samtools_sort_index'
include { MARK_DUPLICATES } from '../../../modules/local/markduplicates'
include { SUBSET_ALIGNMENT } from '../../../modules/local/subset_alignment'

workflow ALT_ALIGN {
    take: 
    ch_fastq
    ref
    fasta_cram

    main:
    BWA_MEM_ALIGN_ALT_POSTALT(
       ref,
       ch_fastq
    )
  
    SAMTOOLS_SORT_INDEX_BEFORE_SUBSET(
                        BWA_MEM_ALIGN_ALT_POSTALT.out.bamfile_postalt,
                        fasta_cram
                        )

    samtools_sorted_index = SAMTOOLS_SORT_INDEX_BEFORE_SUBSET.out.sortedAln

    MARK_DUPLICATES(
        samtools_sorted_index
    )

    SUBSET_ALIGNMENT(
        MARK_DUPLICATES.out.markdupbam,
        fasta_cram
    )

    SAMTOOLS_SORT_INDEX_AFTER_SUBSET(
                                     SUBSET_ALIGNMENT.out.subset_bam,
                                     fasta_cram
                                    )

    emit:
    SAMTOOLS_SORT_INDEX_AFTER_SUBSET.out.sortedAln
}
