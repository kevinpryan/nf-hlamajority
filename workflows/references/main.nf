#!/usr/bin/env nextflow

include { BWA_INDEX as BWA_INDEX_BWAKIT } from '../../modules/local/bwa_index'
include { BWA_INDEX as BWA_INDEX_HLA_LA } from '../../modules/local/bwa_index'
include { BWA_INDEX as BWA_INDEX_KOURAMI } from '../../modules/local/bwa_index'
include { BWA_INDEX as BWA_INDEX_POLYSOLVER } from '../../modules/local/bwa_index'


process GET_IMGT {
    publishDir "${params.references_basedir}/source", mode: 'copy'
    
    label 'HLALA_CONTAINER'
    input:
    val imgt_commit
    output:
    path "IMGTHLA", emit: repo
    
    script:
    """
    git clone https://github.com/ANHIG/IMGTHLA.git
    cd IMGTHLA
    git checkout ${imgt_commit}
    cp ./wmda/hla_nom_g.txt alignments/
    """
}

process GET_KOURAMI {
    publishDir "${params.references_basedir}", mode: 'copy'
    label 'HLALA_CONTAINER'

    input:
    // 545c770
    val kourami_commit
    output:
    path "kourami", emit: repo

    script:
    """
    git clone https://github.com/Kingsford-Group/kourami.git
    cd kourami
    git checkout ${kourami_commit}
    """
}

process KOURAMI_DOWNLOAD_HS38NOALTDH {
    //publishDir "${params.references_basedir}/kourami/resources", mode: 'copy'
    publishDir "${params.references_basedir}/kourami/resources",
               mode: 'copy',
               saveAs: { file -> new File(file).getName() }

    label 'HLALA_CONTAINER'

    input:
    path kourami_repo

    output:
    path("${kourami_repo}/resources/hs38NoAltDH.fa"), emit: reference
    //path("${kourami_repo}/reso/hs38NoAltDH.fa"), emit: reference

    script:
    """
    cd ${kourami_repo}/scripts
    bash download_grch38.sh hs38NoAltDH
    """
}

process BUILD_KOURAMI {
    publishDir "${params.references_basedir}", mode: 'copy'

    container 'kevinr9525/cancerit-kourami:wget' // Nextflow handles the engine switching
    
    input:
    path imgt_repo
    val imgt_version
    path kourami_repo, stageAs: 'kourami_src' 

    output:
    path "kourami" // This captures the directory we create in the script

    script:
    """
    cp -rL ${imgt_repo} imgt_work
    cp -rL kourami_src kourami
    JAR_PATH=/opt/wtsi-cgp/java/Kourami.jar
    cd kourami; mkdir -p target
    cp \$JAR_PATH target
    cd ../
    cp imgt_work/alignments/V_gen.txt imgt_work/alignments/Y_gen.txt
    cp imgt_work/alignments/V_nuc.txt imgt_work/alignments/Y_nuc.txt
    sed -i 's/V\\*/Y\\*/g' imgt_work/alignments/Y_*.txt
    sed -i 's/HLA-V/HLA-Y/g' imgt_work/alignments/Y_*.txt
    cd kourami; mkdir -p custom_db
    cd scripts
    cat formatIMGT.sh > fixed_script.sh
    sed -i 's/Kourami.jar"\\\\];then/Kourami.jar" ];then/' fixed_script.sh
    sed -i 's|cp \$resource_dir/DRB5_gen.txt| echo "removed else condition " # |' fixed_script.sh
    chmod +x fixed_script.sh
    ./fixed_script.sh -i ../../imgt_work/alignments -v ${imgt_version} -o ../custom_db
    """
}

process BUILD_BWAKIT {
    label 'bwa_mem_container'
    publishDir "${params.references_basedir}/bwakit", mode: 'copy'

    output:
    path ("hs38DH.fa"), emit: reference
    path ("hs38DH.fa.alt")

    script:
    """
    run-gen-ref hs38DH
    """
}


process HLA_LA_REFERENCE_DOWNLOAD {
    label 'HLALA_CONTAINER'

    output:
    path("PRG_MHC_GRCh38_withIMGT.tar.gz"), emit: reference_zip

    script:
    """
    wget -O PRG_MHC_GRCh38_withIMGT.tar.gz https://zenodo.org/records/19336310/files/PRG_MHC_GRCh38_withIMGT.tar.gz?download=1 
    """
}

//unstable URL: http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz

process HLA_LA_REFERENCE_PREPARE {
    label 'HLALA_CONTAINER'
    publishDir "${params.references_basedir}/hla-la", mode: 'copy'

    input:
    path(hla_la_reference_gzip)

    output:
    path("PRG_MHC_GRCh38_withIMGT")
    path("PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.fa"), emit: extended_ref

    script:
    """
    tar -xvzf ${hla_la_reference_gzip} && HLA-LA --action prepareGraph --PRG_graph_dir PRG_MHC_GRCh38_withIMGT
    """
}

process POLYSOLVER_REFERENCE_DOWNLOAD {
    label 'HLALA_CONTAINER'
    publishDir "${params.references_basedir}/polysolver", mode: 'copy'

    output:
    path("GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"), emit: reference

    script:
    """
    wget -O "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
    """
}

workflow REFERENCES {
    take:
    reference_dir
    imgt_commit
    imgt_version
    kourami_commit

    main:
    
    GET_IMGT(
            imgt_commit
            )

    GET_KOURAMI(
                kourami_commit
               )

    BUILD_KOURAMI(
                  GET_IMGT.out.repo,
                  imgt_version,
                  GET_KOURAMI.out.repo
                 )

    KOURAMI_DOWNLOAD_HS38NOALTDH(
                                 GET_KOURAMI.out.repo
                                )

    BWA_INDEX_KOURAMI(
                     KOURAMI_DOWNLOAD_HS38NOALTDH.out.reference,
                    "kourami/resources"
                     )

    BUILD_BWAKIT()
  
    BWA_INDEX_BWAKIT(
                     BUILD_BWAKIT.out.reference,
                     "bwakit"
                    )
    Channel
    .from(
        params.hla_la_prg_tar 
            ? file(params.hla_la_prg_tar)
            : null
    )
    .set { ch_hla_la_tar }

    if (params.hla_la_prg_tar) {
        hla_la_zip = file(params.hla_la_prg_tar) 
    } else {
        log.info "No --hla_la_prg_tar provided; performing automated download"
        HLA_LA_REFERENCE_DOWNLOAD()
        hla_la_zip = HLA_LA_REFERENCE_DOWNLOAD.out.reference_zip
    }

    HLA_LA_REFERENCE_PREPARE(hla_la_zip)

    BWA_INDEX_HLA_LA(
                    HLA_LA_REFERENCE_PREPARE.out.extended_ref,
                    "hla-la/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome"
                    )
    
    POLYSOLVER_REFERENCE_DOWNLOAD()

    BWA_INDEX_POLYSOLVER(
                        POLYSOLVER_REFERENCE_DOWNLOAD.out.reference,
                        "polysolver"
                        )

}
