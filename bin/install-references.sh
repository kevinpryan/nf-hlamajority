#!/bin/bash

# Usage: 
#   bash install_references.sh -e docker -o /path/to/refs
#   bash install_references.sh -e singularity -c /path/to/singularity_cache -o /path/to/refs
#
# Arguments:
#   -e, --engine   [Required] Container engine: 'docker' or 'singularity'
#   -c, --cache    [Required if engine is singularity] Directory for singularity images
#   -o, --output   [Optional] Directory to build references (default: ../references)

set -e

# Initialize variables
ENGINE_MODE=""
CACHE_DIR=""
OUTPUT_DIR=""

# ------------------------------------------------------------------
# 1. PARSE NAMED ARGUMENTS
# ------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case $1 in
    -e|--engine)
      ENGINE_MODE="$2"
      shift 2
      ;;
    -c|--cache)
      CACHE_DIR="$2"
      shift 2
      ;;
    -o|--output)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 -e <docker|singularity> [-c <cache_dir>] [-o <output_dir>]"
      echo ""
      echo "Options:"
      echo "  -e, --engine    Container engine (docker or singularity)"
      echo "  -c, --cache     Path to singularity image cache (Required for singularity)"
      echo "  -o, --output    Path to reference output directory (Default: project/references)"
      exit 0
      ;;
    *)
      echo "Error: Unknown argument '$1'"
      echo "Use -h for help."
      exit 1
      ;;
  esac
done

# ------------------------------------------------------------------
# 2. VALIDATE ARGUMENTS
# ------------------------------------------------------------------

# Validate Engine
if [[ -z "$ENGINE_MODE" ]]; then
    echo "Error: You must specify a container engine using -e or --engine."
    echo "       Available options: 'docker', 'singularity'"
    exit 1
fi

if [[ "$ENGINE_MODE" != "docker" && "$ENGINE_MODE" != "singularity" ]]; then
    echo "Error: Invalid engine '$ENGINE_MODE'. Must be 'docker' or 'singularity'."
    exit 1
fi

# Validate Cache (Required for Singularity)
if [[ "$ENGINE_MODE" == "singularity" ]]; then
    if [[ -z "$CACHE_DIR" ]]; then
        echo "Error: When using Singularity, you must specify a cache directory using -c or --cache."
        exit 1
    fi
    # Create and resolve cache path
    mkdir -p "$CACHE_DIR"
    SINGULARITY_CACHE=$(realpath "$CACHE_DIR")
fi

# Validate Output Directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PROJECT_ROOT=$(realpath "${SCRIPT_DIR}/..") 
ASSETS_DIR="${PROJECT_ROOT}/assets"

if [[ -n "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
    REF_DIR=$(realpath "$OUTPUT_DIR")
    echo "INFO: Output directory set to: $REF_DIR"
else
    REF_DIR="${PROJECT_ROOT}/references"
    echo "INFO: Using default output directory: $REF_DIR"
fi

echo "INFO: Engine selected: $ENGINE_MODE"

# ------------------------------------------------------------------
# 3. DEFINE IMAGES
# ------------------------------------------------------------------
# Image IDs: 0=Kourami, 1=BWAKit, 2=HLA-LA
IMAGE_FILES=(
    "kevinr9525-cancerit-kourami-wget.img"
    "quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img"
    "kevinr9525-genome-seek_hla-v1.0.4.img"
)

IMAGE_URIS=(
    "docker://kevinr9525/cancerit-kourami:wget"
    "docker://quay.io/biocontainers/bwakit:0.7.18.dev1--hdfd78af_0"
    "docker://kevinr9525/genome-seek_hla:v1.0.4"
)

# ------------------------------------------------------------------
# 4. PREPARE IMAGES
# ------------------------------------------------------------------
if [ "$ENGINE_MODE" == "singularity" ]; then
    echo "INFO: Verifying Singularity images in cache..."
    for i in "${!IMAGE_FILES[@]}"; do
        FILE="${IMAGE_FILES[$i]}"
        URI="${IMAGE_URIS[$i]}"
        if [ ! -f "$SINGULARITY_CACHE/$FILE" ]; then
            echo "Pulling '$FILE' from '$URI'..."
            singularity pull --name "$SINGULARITY_CACHE/$FILE" "$URI"
        fi
    done
elif [ "$ENGINE_MODE" == "docker" ]; then
    echo "INFO: Ensuring Docker images are pulled..."
    for i in "${!IMAGE_URIS[@]}"; do
        DOCKER_TAG=${IMAGE_URIS[$i]#docker://} 
        # Check if image exists locally
        if [[ "$(docker images -q $DOCKER_TAG 2> /dev/null)" == "" ]]; then
            echo "Pulling $DOCKER_TAG..."
            docker pull "$DOCKER_TAG"
        fi
    done
fi

# ------------------------------------------------------------------
# HELPER: GET CMD
# ------------------------------------------------------------------
get_container_cmd() {
    local idx=$1
    if [ "$ENGINE_MODE" == "singularity" ]; then
        echo "singularity run -B $(pwd) ${SINGULARITY_CACHE}/${IMAGE_FILES[$idx]}"
    else
        local DOCKER_TAG=${IMAGE_URIS[$idx]#docker://}
        # -v $(pwd):$(pwd) maps current dir
        # -w $(pwd) sets work dir
        # -u $(id -u):$(id -g) ensures files are owned by user, not root
        echo "docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) ${DOCKER_TAG}"
    fi
}

# ------------------------------------------------------------------
# HELPER: VERIFY OR RUN
# ------------------------------------------------------------------
verify_or_run() {
    local checksum_file="$1"
    local run_command="$2"
    local step_name="$3"

    echo "============================================================"
    echo "STEP: $step_name"
    
    if [ -f "$checksum_file" ] && md5sum --status -c "$checksum_file"; then
        echo "INFO: Checksum verified. Skipping."
        return 0
    fi

    echo "INFO: Checksum mismatch/missing. Running generation..."
    eval "$run_command"

    if [ -f "$checksum_file" ] && md5sum --status -c "$checksum_file"; then
        echo "SUCCESS: Verified."
    else
        echo "ERROR: Checksum verification failed for $step_name."
        exit 1
    fi
    echo "============================================================"
}

#--------------------------------------------------------------------
# HELPER: TEST HLA-LA INSTALLATION
#--------------------------------------------------------------------

test_hlala() {
    local checksum_file="$1"
    local run_command="$2"
    local step_name="$3"
    
    echo "============================================================"
    echo "STEP: $step_name"
    wget -O NA12878.mini.cram https://www.dropbox.com/scl/fi/7kbq1mrk8k468dscsg3b8/NA12878.mini.cram?rlkey=ng89b7mf913yrk7ppvxcyk20j&e=1&dl=0 
    eval "$run_command"
    # check output of HLA-LA test
    if [ -f "$checksum_file" ] && md5sum --status -c "$checksum_file"; then
       echo "HLA-LA installation successul"
    else 
        echo "ERROR: Checksum verification failed for $step_name - problem with HLA-LA installation or testfile"
        exit 1
    fi
    
}

# ------------------------------------------------------------------
# MAIN INSTALLATION LOGIC
# ------------------------------------------------------------------
mkdir -p ${REF_DIR}; cd ${REF_DIR}

# --- KOURAMI ---
if [ ! -d "kourami" ] ; then
    git clone https://github.com/Kingsford-Group/kourami.git
else
    cd kourami; git pull https://github.com/Kingsford-Group/kourami.git; cd ..
fi

# Kourami: GRCh38

cd "${REF_DIR}/kourami/resources"

verify_or_run \
    "${ASSETS_DIR}/references_checksums/kourami/resources/hs38NoAltDH.fa.md5" \
    "cd ../scripts && bash download_grch38.sh hs38NoAltDH && cd ../resources" \
    "Kourami: Download hs38NoAltDH"

CMD=$(get_container_cmd 0)
verify_or_run \
    "${ASSETS_DIR}/references_checksums/kourami/resources/bwa_0.7.17-r1188_index_hs38NoAltDH.md5" \
    "$CMD bwa index hs38NoAltDH.fa" \
    "Kourami: Index hs38NoAltDH"
cd ../ 

# Kourami: Panel

cd "${REF_DIR}/kourami/"

CMD=$(get_container_cmd 0)
if [ -d "db" ] ; then
    cd db
    verify_or_run \
        "${ASSETS_DIR}/references_checksums/kourami/db/All_FINAL_with_Decoy.fa.gz.md5" \
        "cd .. && $CMD bash scripts/download_panel.sh && cd db" \
        "Kourami: Download Panel"
else
    # need to have db dir to run verify_or_run
    mkdir -p db; cd db
    verify_or_run \
        "${ASSETS_DIR}/references_checksums/kourami/db/All_FINAL_with_Decoy.fa.gz.md5" \
        "cd .. && rm -f db $CMD bash scripts/download_panel.sh && cd db" \
        "Kourami: Download Panel"
fi

verify_or_run \
    "${ASSETS_DIR}/references_checksums/kourami/db/bwa_0.7.18_index_All_Final_with_Decoy.md5" \
    "rm -f All_FINAL_with_Decoy.fa.gz.* && $CMD bwa index All_FINAL_with_Decoy.fa.gz" \
    "Kourami: Index Panel"


# --- BWAKIT ---
cd ${REF_DIR}
mkdir -p bwakit; cd bwakit
CMD=$(get_container_cmd 1)

verify_or_run \
    "${ASSETS_DIR}/references_checksums/bwakit/hs38DH.fa.md5" \
    "rm -f hs38DH.fa hs38DH.fa.alt && $CMD run-gen-ref hs38DH" \
    "BWAKit: Generate Reference"

if ! md5sum --status -c "${ASSETS_DIR}/references_checksums/bwakit/hs38DH.fa.alt.md5"; then
    echo "Error: hs38DH.fa.alt MD5 mismatch."
    exit 1
fi

verify_or_run \
    "${ASSETS_DIR}/references_checksums/bwakit/bwa_0.7.18_index_hs38DH.md5" \
    "$CMD bwa index hs38DH.fa" \
    "BWAKit: Index"


# --- HLA-LA ---
cd ${REF_DIR}
mkdir -p hla-la; cd hla-la
CMD=$(get_container_cmd 2)

verify_or_run \
    "${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT.tar.gz.md5" \
    "rm -f PRG_MHC_GRCh38_withIMGT.tar.gz && wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz" \
    "HLA-LA: Download"

echo "INFO: Building HLA-LA Graph..."

cd ${REF_DIR}/hla-la/ && tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz && $CMD HLA-LA --action prepareGraph --PRG_graph_dir PRG_MHC_GRCh38_withIMGT

echo "INFO: HLA-LA graph built"

# Extended Ref
cd "${REF_DIR}/hla-la/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome"
if [ ! -f "extendedReferenceGenome.fa" ]; then
    echo "Error: extendedReferenceGenome.fa missing."
    exit 1
fi

verify_or_run \
    "${ASSETS_DIR}/references_checksums/hla-la/PRG_MHC_GRCh38_withIMGT/extendedReferenceGenome/bwa_0.7.17-r1188_index_extendedReferenceGenome.md5" \
    "$CMD bwa index extendedReferenceGenome.fa" \
    "HLA-LA: Index Extended Ref"

echo "References completed successfully."
