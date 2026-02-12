include { RUN_HLALA } from '../../../modules/local/run_hlala'

workflow HLA_LA{
    take: 
    bam
    graphdir

    main:
    RUN_HLALA(
        bam,
        graphdir
    )
    emit:
    RUN_HLALA.out.hlala_call
}
