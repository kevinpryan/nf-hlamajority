//
// Alt-aware alignment with bwa-mem and postprocessing script
// Also subsets output bam file to relevant regions for HLA typing
//

//include {bwa_mem_align_alt_postalt} from '../../../modules/local/bwa_mem_align_alt_postalt'
include {samtools_sort} from '../../../modules/local/samtools_sort'
include {samtools_index} from '../../../modules/local/samtools_index'
include {markduplicates} from '../../../modules/local/markduplicates'
include {extractContigs} from '../../../modules/local/extractContigs'
include {fasta_index_bed} from '../../../modules/local/fasta_index_bed'
include {subsetBam} from '../../../modules/local/subsetBam'

workflow align{
    take: 
    ch_fastq
    ch_ref
    ch_hlatypes
    reference_basename
    chr
    main:
    bwa_mem_align(
       ch_ref,
       ch_fastq,
       reference_basename
    )  
    samtools_sort(
        bwa_mem_align_alt_postalt.out.bamfile
    )
    samtools_index(
        samtools_sort.out.sortedbam
    )
    samtools_sorted_index = samtools_sort.out.sortedbam.join(samtools_index.out.bam_indexed, by: 0)
    markduplicates(
        samtools_sorted_index
    )
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
        fasta_index_bed.out.fasta_bed
    )
    emit: 
    subsetBam.out.subsetbam
}

