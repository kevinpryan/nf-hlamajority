include { RUN_HLALA } from '../../../modules/local/run_hlala'
include { RUN_HLALA_PLACEHOLDER } from '../../../modules/local/run_hlala'

workflow HLA_LA {
    take: 
    bam
    graphdir

    main:
    RUN_HLALA(
        bam,
        graphdir
    )
    // need to deal with two scenarios: 1. the user provides a single-end file (skip RUN_HLALA_PLACEHOLDER) 2. RUN_HLALA fails (ignore failure)
    // Determine who succeeded and who failed/skipped
    bam
        // Create a key-only channel [meta] to match against results
        //.map { meta, bam, bai -> meta }
        .map { meta, bam, bai -> [ meta ] } // wrap in list to make it a tuple key
        
        // Join with results. 'remainder: true' keeps keys even if results are missing.
        // Result format: [ meta, hlala_output_or_null ]
        .join(RUN_HLALA.out.hlala_call, remainder: true)
        
        // Split into two channels based on whether we got a result
        .branch { meta, results ->
            success: results != null
                return [meta, results]  // Pass the file path
            failure: results == null
                return meta    // Pass just the meta for the placeholder
        }
        .set { ch_routing }

    // Run placeholder for EVERYONE who didn't get a result
    // This includes Single-End (skipped) AND Paired-End (failed/ignored)
    RUN_HLALA_PLACEHOLDER(
        ch_routing.failure
    )

    // Merge results
    // We reconstruct the [meta, path] tuple for the success channel 
    // to match the placeholder output structure if needed, 
    // or just mix if structures align.
    
    // ch_routing.success looks like [path] or [meta, path] 
    // depending on exact join behavior. 
    // join(remainder:true) returns [key, val]. 
    // So ch_routing.success is [meta, path].

    emit:
    calls = ch_routing.success.mix(RUN_HLALA_PLACEHOLDER.out.hlala_call)
}
