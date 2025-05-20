include { RUN_HLALA } from '../../../modules/local/run_hlala'

workflow hlala{
    /*
    convert bam to fastq then run Optitype
    */
    take: 
    subsetbam
    graphdir

    main:
    RUN_HLALA(
        subsetbam,
        graphdir
    )
    emit:
    RUN_HLALA.out.hlala_call
}
