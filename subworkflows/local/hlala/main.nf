include { RUN_HLALA } from '../../../modules/local/run_hlala'
include { RUN_HLALA_PLACEHOLDER_SINGLE_END } from '../../../modules/local/run_hlala'
include { RUN_HLALA_PLACEHOLDER_FAILURE } from '../../../modules/local/run_hlala'

workflow HLA_LA {

    take:
    bam
    graphdir

    main:

    RUN_HLALA(
        bam,
        graphdir
    )


    // Expected samples: keep only key
    bam
        .map { meta, bam, index ->
            [meta.sample, meta]
        }
        .set { expected_samples }


    // Successful HLA-LA outputs
    RUN_HLALA.out.hlala_call
        .map { meta, result ->
            [meta.sample, result]
        }
        .set { successful_hlala }


    // Find missing samples
    expected_samples
        .join(successful_hlala, remainder:true)
        .filter { sample, meta, result ->
            result == null
        }
        .map { sample, meta, result ->
            meta
        }
        .branch {
            paired: it.single_end == false
            single: it.single_end == true
        }
        .set { failed_hlala }

    failed_hlala_single = failed_hlala.single
    failed_hlala_paired = failed_hlala.paired

    RUN_HLALA_PLACEHOLDER_SINGLE_END(
        failed_hlala.single
    )

    RUN_HLALA_PLACEHOLDER_FAILURE(
        failed_hlala.paired
    )
    
    emit:

    hlala_call =
        RUN_HLALA.out.hlala_call
        .mix(RUN_HLALA_PLACEHOLDER_SINGLE_END.out.hlala_call)
        .mix(RUN_HLALA_PLACEHOLDER_FAILURE.out.hlala_call)
   
    status = RUN_HLALA.out.run_status
            .mix(RUN_HLALA_PLACEHOLDER_SINGLE_END.out.run_status)
            .mix(RUN_HLALA_PLACEHOLDER_FAILURE.out.run_status)
    
}
