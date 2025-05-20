include { bam2fastq } from '../../../modules/local/bam2fastq'
include { RUN_OPTITYPE } from '../../../modules/local/run_optitype'

workflow optitype{
    /*
    convert bam to fastq then run Optitype
    */
    take: 
    bam
    //dna_rna

    main:
    bam2fastq(
        bam
    )
    RUN_OPTITYPE(
        bam2fastq.out.convertedfastqs//,
        //dna_rna
    )
    emit:
    RUN_OPTITYPE.out.optitype_call
}
