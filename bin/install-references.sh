#!/bin/bash

# to run: go to bin directory
# bash install_references.sh /path/to/dir/with/singularity/images
# HLA-LA prepareGraph can take a few hours and take up to 40G of memory
set -e

if [ -z "$1" ]; then
    echo "Error: Please provide the path to the singularity cache directory as the first argument."
    echo "Usage: $0 <singularity_cache_dir> [reference_output_dir]"
    echo "Example: $0 /path/to/singularity_images /data/large_storage/references"
    exit 1
fi

SINGULARITY_CACHE=$1

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PROJECT_ROOT=$(realpath "${SCRIPT_DIR}/..") # Assumes script is in project/bin
ASSETS_DIR="${PROJECT_ROOT}/assets"

if [ -n "$2" ]; then
    # User provided a custom path
    # Create it immediately to ensure realpath works and permissions are okay
    mkdir -p "$2"
    REF_DIR=$(realpath "$2")
    echo "INFO: Custom reference directory specified: $REF_DIR"
else
    # Default to project/references
    REF_DIR="${PROJECT_ROOT}/references"
    echo "INFO: Using default reference directory: $REF_DIR"
fi

# ------------------------------------------------------------------
# HELPER FUNCTION: Verify MD5, Run if missing/bad, Verify again
# ------------------------------------------------------------------
verify_or_run() {
    local checksum_file="$1"
    local run_command="$2"
    local step_name="$3"

    echo "============================================================"
    echo "STEP: $step_name"
    
    # Check 1: Do files exist and match checksum?
    # We suppress output here to avoid clutter if it passes
    if [ -f "$checksum_file" ] && md5sum --status -c "$checksum_file"; then
        echo "INFO: Checksum verified. Skipping installation for $step_name."
        return 0
    fi

    echo "INFO: Checksum missing or mismatch. Running installation..."
    
    # Run the command
    eval "$run_command"

    # Check 2: Verify after installation
    if [ -f "$checksum_file" ] && md5sum --status -c "$checksum_file"; then
        echo "SUCCESS: $step_name installed and verified successfully."
    else
        echo "ERROR: $step_name failed verification after installation."
        echo "       Checked against: $checksum_file"
        exit 1
    fi
    echo "============================================================"
}
# ------------------------------------------------------------------


mkdir -p ${SINGULARITY_CACHE}; cd $SINGULARITY_CACHE

# Pull singularity images if they are not already present
IMAGE_FILES=(
    "kevinr9525-cancerit-kourami-wget.img"
    "quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img"
    "kevinr9525-genome-seek_hla-v1.0.4.img"
)

# Define the corresponding remote Docker URIs
IMAGE_URIS=(
    "docker://kevinr9525/cancerit-kourami:wget"
    "docker://quay.io/biocontainers/bwakit:0.7.18.dev1--hdfd78af_0"
    "docker://kevinr9525/genome-seek_hla:v1.0.4"
)

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

# ------------------------------------------------------------------
# 1. BUILD/DOWNLOAD KOURAMI REFERENCES
# ------------------------------------------------------------------
mkdir -p ${REF_DIR}; cd ${REF_DIR}

# Ensure repo is present
if [ ! -d "kourami" ] ; then
    git clone https://github.com/Kingsford-Group/kourami.git
else
    # Verify git repo is valid, otherwise re-clone could be better, but pull is standard
    cd kourami; git pull https://github.com/Kingsford-Group/kourami.git; cd ..
fi

# A. Download GRCh38
# We verify inside 'resources' because md5sum checks relative paths
cd kourami/resources
verify_or_run \
    "${ASSETS_DIR}/references_checksums/kourami/resources/hs38NoAltDH.fa.md5" \
    "cd ../scripts && bash download_grch38.sh hs38NoAltDH && cd ../resources" \
    "Kourami: Download hs38NoAltDH"

# B. Index GRCh38
verify_or_run \
    "${ASSETS_DIR}/references_checksums/kourami/resources/bwa_0.7.17-r1188_index_hs38NoAltDH.md5" \
    "singularity run -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-cancerit-kourami-wget.img bwa index hs38NoAltDH.fa" \
    "Kourami: Index hs38NoAltDH"

cd ../ # back to kourami root

# C. Download Panel
# We verify inside 'db'
cd db
verify_or_run \
    "${ASSETS_DIR}/references_checksums/kourami/db/All_FINAL_with_Decoy.fa.gz.md5" \
    "cd .. && singularity run --cleanenv -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-cancerit-kourami-wget.img bash scripts/download_panel.sh && cd db" \
    "Kourami: Download Panel"

# D. Index Panel
verify_or_run \
    "${ASSETS_DIR}/references_checksums/kourami/db/bwa_0.7.18_index_All_Final_with_Decoy.md5" \
    "rm -f All_FINAL_with_Decoy.fa.gz.* && singularity run -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-cancerit-kourami-wget.img bwa index All_FINAL_with_Decoy.fa.gz" \
    "Kourami: Index Panel"


# ------------------------------------------------------------------
# 2. BUILD/DOWNLOAD BWAKIT REFERENCES
# ------------------------------------------------------------------
cd ${REF_DIR}
mkdir -p bwakit; cd bwakit

# A. Run Gen Ref
# Note: Checking both fa and alt md5s
verify_or_run \
    "${ASSETS_DIR}/references_checksums/bwakit/hs38DH.fa.md5" \
    "rm -f hs38DH.fa hs38DH.fa.alt && singularity run -B $(pwd) ${SINGULARITY_CACHE}/quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img run-gen-ref hs38DH" \
    "BWAKit: Generate Reference (hs38DH)"

# Double check .alt file specifically just in case
if ! md5sum --status -c "${ASSETS_DIR}/references_checksums/bwakit/hs38DH.fa.alt.md5"; then
    echo "Error: hs38DH.fa.alt MD5 mismatch even after generation."
    exit 1
fi

# B. Index BWAKit
verify_or_run \
    "${ASSETS_DIR}/references_checksums/bwakit/bwa_0.7.18_index_hs38DH.md5" \
    "singularity run -B $(pwd) ${SINGULARITY_CACHE}/quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img bwa index hs38DH.fa" \
    "BWAKit: Index hs38DH"


# ------------------------------------------------------------------
# 3. BUILD/DOWNLOAD HLA-LA REFERENCES
# ------------------------------------------------------------------
cd ${REF_DIR}
mkdir -p hla-la; cd hla-la

# A. Download Tarball
verify_or_run \
    "${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT.tar.gz.md5" \
    "rm -f PRG_MHC_GRCh38_withIMGT.tar.gz && wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz" \
    "HLA-LA: Download Data Package"

# B. Extract and Prepare Graph
# We check the Final Output (serializedGRAPH) to see if we can skip the extraction/prep
cd PRG_MHC_GRCh38_withIMGT
# Check if graph exists
if [ -f "${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/serializedGRAPH.md5" ] && \
   md5sum --status -c "${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/serializedGRAPH.md5"; then
   
   echo "INFO: HLA-LA Graph already prepared and verified. Skipping extraction and graph prep."

else
   # If graph is missing/bad, we must go back up, extract, and run
   echo "INFO: HLA-LA Graph missing or invalid. Extracting and building..."
   cd .. # Back to hla-la dir where tar.gz is
   
   # Untar
   tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz
   # We DO NOT remove the tarball yet, in case the next step fails and we need to retry later
   
   # Prepare Graph
   singularity run -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-genome-seek_hla-v1.0.4.img HLA-LA --action prepareGraph --PRG_graph_dir PRG_MHC_GRCh38_withIMGT
   
   # Now we verify the result
   cd PRG_MHC_GRCh38_withIMGT
   if md5sum --status -c "${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/serializedGRAPH.md5"; then
       echo "SUCCESS: HLA-LA Graph prepared successfully."
       # Safe to clean up tarball now
       rm ../PRG_MHC_GRCh38_withIMGT.tar.gz
   else 
       echo "ERROR: HLA-LA Graph generation failed checksum."
       exit 1
   fi
fi

# C. Index Extended Reference Genome
cd extendedReferenceGenome
verify_or_run \
    "${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/extendedReferenceGenome.fa.md5" \
    "echo 'Error: Extended reference genome FASTA should have been extracted by previous steps. Check Graph Prep.' && exit 1" \
    "HLA-LA: Check Extended Reference Fasta"

verify_or_run \
    "${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/bwa_0.7.17-r1188_index_extendedReferenceGenome.md5" \
    "singularity run -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-genome-seek_hla-v1.0.4.img bwa index extendedReferenceGenome.fa" \
    "HLA-LA: Index Extended Reference"

echo "References have been downloaded and built."
