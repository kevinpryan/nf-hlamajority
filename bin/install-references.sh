#!/bin/bash

# to run: go to bin directory
# bash install_references.sh /path/to/dir/with/singularity/images
# HLA-LA prepareGraph can take a few hours and take up to 40G of memory
set -e

if [ -z "$1" ]; then
    echo "Error: Please provide the path to the singularity cache directory as the first argument."
    echo "Usage: $0 /path/to/singularity/images"
    exit 1
fi

SINGULARITY_CACHE=$1

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PROJECT_ROOT=$(realpath "${SCRIPT_DIR}/..") # Assumes script is in project/bin
REF_DIR="${PROJECT_ROOT}/references"
ASSETS_DIR="${PROJECT_ROOT}/assets"

#bindir=$(pwd)

mkdir -p ${SINGULARITY_CACHE}; cd $SINGULARITY_CACHE

# Pull singularity images if they are not already present
IMAGE_FILES=(
    "kevinr9525-cancerit-kourami-wget.img"
    "quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img"
    "kevinr9525-genome-seek_hla-v1.0.4.img"
)
#     "bwa%3A0.7.18--he4a0461_1"

# Define the corresponding remote Docker URIs
IMAGE_URIS=(
    "docker://kevinr9525/cancerit-kourami:wget"
    "docker://quay.io/biocontainers/bwakit:0.7.18.dev1--hdfd78af_0"
    "docker://kevinr9525/genome-seek_hla:v1.0.4"
)
#     "https://depot.galaxyproject.org/singularity/bwa%3A0.7.18--he4a0461_1"
# Loop over the images by their index
for i in "${!IMAGE_FILES[@]}"; do
    FILE="${IMAGE_FILES[$i]}"
    URI="${IMAGE_URIS[$i]}"
    if [ ! -f "$FILE" ]; then
        echo "Image not found. Pulling '$FILE' from '$URI'..."
        singularity pull --name "$FILE" "$URI"
    else
        echo "Image '$FILE' already exists. Using existing image."
    fi
done

# BUILD/DOWNLOAD KOURAMI REFERENCES
mkdir -p ${REF_DIR}; cd ${REF_DIR}
git clone https://github.com/Kingsford-Group/kourami.git
cd kourami/scripts
bash download_grch38.sh hs38NoAltDH
cd ../resources/
if [ -f ${ASSETS_DIR}/references_checksums/kourami/resources/hs38NoAltDH.fa.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/kourami/resources/hs38NoAltDH.fa.md5; then
    echo "INFO: Reference genome correctly downloaded. Indexing."
    singularity run -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-cancerit-kourami-wget.img bwa index hs38NoAltDH.fa
else
    echo "INFO: Reference genome kourami not correctly downloaded. Exiting" >&2
    exit 1
fi

if [ -f ${ASSETS_DIR}/references_checksums/kourami/resources/bwa_0.7.17-r1188_index_hs38NoAltDH.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/kourami/resources/bwa_0.7.17-r1188_index_hs38NoAltDH.md5; then
    echo "Kourami reference indexing successful"
else
    echo "Kourami reference indexing unsuccessful. Exiting." >&2
    exit 1
fi

cd ../

# cannot use kourami image for this as it doesn't contain
#singularity run --cleanenv -B $(pwd) ${SINGULARITY_CACHE}/bwa%3A0.7.18--he4a0461_1 bash scripts/download_panel.sh
# kourami image contains bwa version 0.7.17-r1188 
singularity run --cleanenv -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-cancerit-kourami-wget.img bash scripts/download_panel.sh
cd db
if [ -f ${ASSETS_DIR}/references_checksums/kourami/db/All_FINAL_with_Decoy.fa.gz.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/kourami/db/All_FINAL_with_Decoy.fa.gz.md5; then
    echo "INFO: Kourami panel downloaded correctly."
else
    echo "INFO: Kourami panel not downloaded correctly. Exiting." >&2
    exit 1
fi

if [ -f ${ASSETS_DIR}/references_checksums/kourami/db/bwa_0.7.18_index_All_Final_with_Decoy.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/kourami/db/bwa_0.7.18_index_All_Final_with_Decoy.md5; then
    echo "INFO: Kourami panel has been indexed correctly."
else
   echo "INFO: Kourami panel not indexed correctly. Exiting." >&2
   exit 1
fi

# BUILD/DOWNLOAD BWAKIT REFERENCES

cd ${REF_DIR}
mkdir -p bwakit; cd bwakit
singularity run -B $(pwd) ${SINGULARITY_CACHE}/quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img run-gen-ref hs38DH
if [ -f ${ASSETS_DIR}/references_checksums/bwakit/hs38DH.fa.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/bwakit/hs38DH.fa.md5 && [ -f ${ASSETS_DIR}/references_checksums/bwakit/hs38DH.fa.alt.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/bwakit/hs38DH.fa.alt.md5; then
    echo "INFO: bwakit references (hs38DH.fa and/or hs38DH.fa.alt) downloaded correctly. Indexing"
else
    echo "INFO: bwakit references (hs38DH.fa and/or hs38DH.fa.alt) not downloaded correctly. Exiting." >&2
    exit 1
fi

singularity run -B $(pwd) ${SINGULARITY_CACHE}/quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img bwa index hs38DH.fa
if  [ -f ${ASSETS_DIR}/references_checksums/bwakit/bwa_0.7.18_index_hs38DH.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/bwakit/bwa_0.7.18_index_hs38DH.md5; then
    echo "INFO: bwakit hs38DH has been indexed correctly."
else
   echo "INFO: bwakit hs38DH has not indexed correctly. Exiting." >&2
   exit 1
fi

# BUILD/DOWNLOAD HLA-LA REFERENCES
cd ${REF_DIR}
mkdir -p hla-la; cd hla-la
wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz

if  [ -f ${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT.tar.gz.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT.tar.gz.md5; then
    echo "INFO: HLA-LA data package has downloaded correctly. Proceeding to unzip and remove PRG_MHC_GRCh38_withIMGT.tar.gz."
else
    echo "INFO: HLA-LA data package has not downloaded correctly. Exiting." >&2
    exit 1
fi

tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz
rm PRG_MHC_GRCh38_withIMGT.tar.gz

singularity run -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-genome-seek_hla-v1.0.4.img HLA-LA --action prepareGraph --PRG_graph_dir PRG_MHC_GRCh38_withIMGT
cd PRG_MHC_GRCh38_withIMGT
if [ -f ${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/serializedGRAPH.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/serializedGRAPH.md5 && [ -f ${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/serializedGRAPH_preGapPathIndex.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/serializedGRAPH_preGapPathIndex.md5; then
    echo "INFO: HLA-LA graph has been prepared correctly"
else
    echo "INFO: HLA-LA graph has not been prepared correctly. Exiting." >&2
    exit 1
fi

#echo "contents of hla-la extendedreferences directory..."
#du -h PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.fa*
cd extendedReferenceGenome/
if  [ -f ${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.fa.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.fa.md5; then
    echo "INFO: HLA-LA extended reference genome has been downloaded correctly. Indexing."
    singularity run -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-genome-seek_hla-v1.0.4.img bwa index extendedReferenceGenome.fa
else
    echo "INFO: HLA-LA extended reference genome has not been downloaded correctly. Exiting." >&2
    exit 1
fi

if  [ -f ${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/bwa_0.7.17-r1188_index_extendedReferenceGenome.md5 ] && md5sum --status -c ${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/bwa_0.7.17-r1188_index_extendedReferenceGenome.md5; then
    echo "INFO: HLA-LA extended reference genome has been indexed correctly."
else
    echo "INFO: HLA-LA extended reference genome has not been index correctly. Exiting." >&2
    exit 1
fi

echo "References have been downloaded and built."
