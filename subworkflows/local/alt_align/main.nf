//
// Alt-aware alignment with bwa-mem and postprocessing script
// Also subsets output bam file to relevant regions for HLA typing
//

include {bwa_mem_align_alt_postalt} from '../../../modules/local/bwa_mem_align_alt_postalt'
//include {samtools_sort} from '../../../modules/local/samtools_sort'
//include {samtools_index} from '../../../modules/local/samtools_index'
include { samtools_sort_index as samtools_sort_index_before_subset } from '../../../modules/local/samtools_sort_index'
include { samtools_sort_index as samtools_sort_index_after_subset } from '../../../modules/local/samtools_sort_index'
include {markduplicates} from '../../../modules/local/markduplicates'
include {extractContigs} from '../../../modules/local/extractContigs'
include {fasta_index_bed} from '../../../modules/local/fasta_index_bed'
//include {subsetBam2} from '../../../modules/local/subsetBam'
include { SUBSET_ALIGNMENT} from '../../../modules/local/subset_alignment'

workflow alt_align{
    take: 
    ch_fastq
    ch_ref
    ch_fasta_cram

    main:
    bwa_mem_align_alt_postalt(
       ch_ref,
       ch_fastq//,
       //reference_basename
    )  
    /*
    samtools_sort(
        bwa_mem_align_alt_postalt.out.bamfile_postalt
    )
    samtools_index(
        samtools_sort.out.sortedbam
    )
    */
    samtools_sort_index_before_subset(
                        bwa_mem_align_alt_postalt.out.bamfile_postalt,
                        ch_fasta_cram
                        )

    //samtools_sorted_index = samtools_sort.out.sortedbam.join(samtools_index.out.bam_indexed, by: 0)
    samtools_sorted_index = samtools_sort_index_before_subset.out.sortedAln

    markduplicates(
        samtools_sorted_index
    )
    /*
    subsetBam2(
        markduplicates.out.markdupbam
    )
    */
    SUBSET_ALIGNMENT(
        markduplicates.out.markdupbam,
        ch_fasta_cram
    )
    samtools_sort_index_after_subset(
                                     SUBSET_ALIGNMENT.out.subset_bam,
                                     ch_fasta_cram
                                     )
                                     

    emit:
    samtools_sort_index_after_subset.out.sortedAln
}

