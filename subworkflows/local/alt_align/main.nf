//
// Alt-aware alignment with bwa-mem and postprocessing script
// Also subsets output bam file to relevant regions for HLA typing
//

include {bwa_mem_align_alt_postalt} from '../../../modules/local/bwa_mem_align_alt_postalt'
include {samtools_sort} from '../../../modules/local/samtools_sort'
include {samtools_index} from '../../../modules/local/samtools_index'
include {markduplicates} from '../../../modules/local/markduplicates'
include {extractContigs} from '../../../modules/local/extractContigs'
include {fasta_index_bed} from '../../../modules/local/fasta_index_bed'
include {subsetBam} from '../../../modules/local/subsetBam'
include {subsetBam2} from '../../../modules/local/subsetBam'


workflow alt_align{
    take: 
    ch_fastq
    ch_ref
    ch_hlatypes
    reference_basename
    chr
    ch_subset_regions
    main:
    bwa_mem_align_alt_postalt(
       ch_ref,
       ch_fastq,
       reference_basename
    )  
    samtools_sort(
        bwa_mem_align_alt_postalt.out.bamfile_postalt
    )
    samtools_index(
        samtools_sort.out.sortedbam
    )
    samtools_sorted_index = samtools_sort.out.sortedbam.join(samtools_index.out.bam_indexed, by: 0)
    markduplicates(
        samtools_sorted_index
    )
// comment out starting here to remove subsetbam
    /*
    extractContigs(
        ch_hlatypes,
        ch_ref,
        reference_basename,
        chr
    )
    fasta_index_bed(
        ch_ref,
        chr,
        reference_basename
    )
    subsetBam(
        markduplicates.out.markdupbam,
        extractContigs.out.alt_contigs,
        extractContigs.out.hla_contigs,
        ch_ref,
        fasta_index_bed.out.fasta_bed,
        ch_subset_regions
    )
*/


   subsetBam2(
        markduplicates.out.markdupbam
    )

// comment out ending here to remove subsetbam
    emit:
    subsetBam2.out.subsetbam 
    //subsetBam.out.subsetbam
    //markduplicates.out.markdupbam
}

