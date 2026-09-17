//
// takes subset bam file from alt_align workflow and aligns to kourami reference. then carry out kourami hla typing
// 

include { RUN_KOURAMI_ALIGN_EXTRACT } from '../../../modules/local/run_kourami/align_extract'
include { RUN_KOURAMI_JAR } from '../../../modules/local/run_kourami/kourami_jar'
include { RUN_KOURAMI_PLACEHOLDER } from '../../../modules/local/run_kourami/kourami_jar'
include { RUN_KOURAMI_PLACEHOLDER_SE } from '../../../modules/local/run_kourami/align_extract'

workflow KOURAMI {
    
    take: 
    bam_ch
    kourami_panel
    kourami_reference

    main:
    RUN_KOURAMI_ALIGN_EXTRACT(
        bam_ch,
        kourami_panel,
        kourami_reference
    )
    
    // Expected samples: keep only key
    bam_ch
        .map { meta, bam, index ->
            [meta.sample, meta]
        }
        .set { expected_samples }
     
    // Successful Kourami align extract outputs
    RUN_KOURAMI_ALIGN_EXTRACT.out.kourami_alignment
        .map { meta, bam_kourami ->
            [meta.sample, bam_kourami]
        }
        .set { successful_kourami_align_extract }

    // Find missing samples
    expected_samples
        .join(successful_kourami_align_extract, remainder:true)
        .filter { sample, meta, bam_kourami ->
            bam_kourami == null
        }
        .map { sample, meta, bam_kourami ->
            meta
        }
        .branch {
            paired: it.single_end == false
            single: it.single_end == true
        }
        .set { failed_kourami_align_extract }
    
    RUN_KOURAMI_PLACEHOLDER_SE(failed_kourami_align_extract.single)

    RUN_KOURAMI_ALIGN_EXTRACT.out.kourami_alignment
        .filter { meta, bam_kourami ->
            bam_kourami != null    
    }
    .set { successful_kourami_align_extract } 

    RUN_KOURAMI_JAR(
        successful_kourami_align_extract,
        kourami_panel
    )

    // Expected samples Kourami Jar
    successful_kourami_align_extract
       .map { meta, bam_kourami ->
            [meta.sample, meta]
       }
       .set { expected_samples_kourami_jar }
    
    // Success Kourami jar extract outputs
    RUN_KOURAMI_JAR.out.kourami_result
        .map { meta, kourami_result ->
            [meta.sample, kourami_result]
        }
        .set { successful_kourami_jar }

    // Find missing samples Kourami jar
    expected_samples_kourami_jar
        .join(successful_kourami_jar, remainder: true)
        .filter{ sample, meta, kourami_result ->
            kourami_result == null
    }
    .map { sample, meta, kourami_result ->
           meta
 
         }
         .set { failed_kourami_jar }

    RUN_KOURAMI_PLACEHOLDER( failed_kourami_jar )

    emit:
    calls = RUN_KOURAMI_JAR.out.kourami_result.mix(RUN_KOURAMI_PLACEHOLDER.out.kourami_result)
                                              .mix(RUN_KOURAMI_PLACEHOLDER_SE.out.kourami_result)
    status = RUN_KOURAMI_JAR.out.run_status
            .mix(RUN_KOURAMI_PLACEHOLDER.out.run_status)
            .mix(RUN_KOURAMI_PLACEHOLDER_SE.out.run_status) 
}
