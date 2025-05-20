//
// takes subset bam file from alt_align workflow and aligns to kourami reference. then carry out kourami hla typing
// 

include { RUN_KOURAMI_ALIGN_EXTRACT } from '../../../modules/local/run_kourami/align_extract'
include { RUN_KOURAMI_JAR } from '../../../modules/local/run_kourami/kourami_jar'

workflow kourami{
    
    take: 
    subsetbam
    kourami_panel
    kourami_reference

    main:
    RUN_KOURAMI_ALIGN_EXTRACT(
        subsetbam,
        kourami_panel,
        kourami_reference
    )

    RUN_KOURAMI_JAR(
        RUN_KOURAMI_ALIGN_EXTRACT.out.kourami_alignment,
        kourami_panel
    )
    emit:
    RUN_KOURAMI_JAR.out.kourami_result
}
