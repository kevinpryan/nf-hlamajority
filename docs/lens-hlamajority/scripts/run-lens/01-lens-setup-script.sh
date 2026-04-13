#!/bin/bash
set -x

usage() {
    echo "Usage: $0 [-p <prefix> -m <manifest> -r <raft_dir>]" 1>&2
    echo "  -p : project name"
    echo "  -m : manifest name (must be on S3 bucket under metadata/ots)"
    echo "  -r : absolute path to directory where raft is installed"
    echo "  -b : base path to S3 bucket where references are located. References are at s3://bucket-base/references/cloud/. nf-hlamajority references are at s3://bucket-base/references/cloud/nf-hlamajority/"
    exit 1
}

while getopts ":p:m:r:b:" o; do
    case "${o}" in
        p) p=${OPTARG} ;;
	m) m=${OPTARG} ;;
	r) r=${OPTARG} ;;
	b) b=${OPTARG} ;;
        *) usage ;;
    esac
done
shift $((OPTIND-1))

if [[ -z "${p}" || -z "${m}" || -z "${r}" || -z "${b}" ]]; then
    echo "something went wrong passing arguments"
    usage
fi

start_time=$SECONDS

RAFT_PATH=${r}
cd $RAFT_PATH

projectDir=${p}
manifest=${m}
export LENS_PROJ_ID=${projectDir}

bucket=${b}

# set up project with manifest
# doesn't run it, just sets up and loads most of the references

###############################
##                           ##
## STEP ONE: SET UP PROJECT  ##
##                           ## 
###############################

raft.py run-ots --project-id $LENS_PROJ_ID \
                --workflow lens \
                --version 1.7-dev-hlamajority \
                --species human \
                --input-data fastqs \
                --manifest $manifest \
                --cloud ${bucket}/references/cloud \
		--setup-only

###############################
##                           ##
## STEP TWO: ALTER CONFIG    ##
##                           ##
###############################

CONFIG_FILE="${RAFT_PATH}/nextflow-for-aws-lens.config"
mkdir -p "${RAFT_PATH}/setup-lens-hlamajority/${projectDir}"
CONFIG_FILE_COMPLETE="${RAFT_PATH}/setup-lens-hlamajority/${projectDir}/nextflow-for-aws-lens-${LENS_PROJ_ID}-complete.config"
cp $CONFIG_FILE $CONFIG_FILE_COMPLETE
mkdir -p projects/${LENS_PROJ_ID}/bams
# overwrite mounts.config with new docker-friendly format
cd ${RAFT_PATH}/projects/${LENS_PROJ_ID}/workflow/
awk -F, '{ for(i=1; i<=NF; i++) { printf "-v %s:%s%s", $i, $i, (i==NF ? "" : " ") } }' mounts.config > mounts.config.tmp && mv mounts.config.tmp mounts.config

# Append the dynamic containerOptions to the AWS config file
cat << EOF >> ${CONFIG_FILE_COMPLETE}

// Appended by launch script to mount project-specific paths
process {
    containerOptions = { "\`cat ${RAFT_PATH}/projects/${LENS_PROJ_ID}/workflow/mounts.config\`" }
}
EOF

echo "Configuration complete. Ready to launch Nextflow."

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time 01-lens-setup.sh: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
