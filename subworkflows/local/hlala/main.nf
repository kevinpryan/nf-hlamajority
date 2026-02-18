include { RUN_HLALA } from '../../../modules/local/run_hlala'
include { RUN_HLALA_PLACEHOLDER } from '../../../modules/local/run_hlala'

workflow HLA_LA{
    take: 
    bam
    graphdir

    main:
    RUN_HLALA(
        bam,
        graphdir
    )
    
    RUN_HLALA_PLACEHOLDER(
        bam.map { meta, aln, bai -> meta }
    )


    emit:
    //RUN_HLALA.out.hlala_call
 RUN_HLALA.out.hlala_call
        .mix(RUN_HLALA_PLACEHOLDER.out.hlala_call)
}
