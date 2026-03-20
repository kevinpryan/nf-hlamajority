//
// takes subset bam file from alt_align workflow and aligns to kourami reference. then carry out kourami hla typing
// 

include { RUN_KOURAMI_ALIGN_EXTRACT } from '../../../modules/local/run_kourami/align_extract'
include { RUN_KOURAMI_JAR } from '../../../modules/local/run_kourami/kourami_jar'
include { RUN_KOURAMI_PLACEHOLDER } from '../../../modules/local/run_kourami/kourami_jar'

workflow KOURAMI {
    
    take: 
    bam
    kourami_panel
    kourami_reference

    main:
    RUN_KOURAMI_ALIGN_EXTRACT(
        bam,
        kourami_panel,
        kourami_reference
    )

    RUN_KOURAMI_JAR(
        RUN_KOURAMI_ALIGN_EXTRACT.out.kourami_alignment,
        kourami_panel
    )

    bam
        // Create a key-only channel from input [meta]
        //.map { meta, bam, bai -> meta }
        .map { meta, bam, bai -> [ meta ] } // wrap in list to make it a tuple key
        // Join with output. If output is missing (timeout), result is null.
        .join(RUN_KOURAMI_JAR.out.kourami_result, remainder: true)
        
        // Split into Success vs Failure
        .branch { meta, result ->
            success: result != null
                return [meta, result] // Return the [meta, path] tuple
            failure: result == null
                return meta   // Return just [meta] for the placeholder
        }
        .set { ch_kourami_routing }

    // 3. Run Placeholder for Timed Out samples
    RUN_KOURAMI_PLACEHOLDER(
        ch_kourami_routing.failure
    )

    emit:
    //RUN_KOURAMI_JAR.out.kourami_result
    calls = ch_kourami_routing.success.mix(RUN_KOURAMI_PLACEHOLDER.out.kourami_result)
}
