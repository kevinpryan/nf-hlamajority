include { BAM_TO_FASTQ } from '../../../modules/local/bam_to_fastq'
include { RUN_OPTITYPE } from '../../../modules/local/run_optitype'

workflow OPTITYPE {
    /*
    convert bam to fastq then run Optitype
    */
    take: 
    bam

    main:
    BAM_TO_FASTQ(
        bam
    )
    RUN_OPTITYPE(
        BAM_TO_FASTQ.out.convertedfastqs
    )

    emit:
    RUN_OPTITYPE.out.optitype_call
}

