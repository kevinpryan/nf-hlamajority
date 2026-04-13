#!/bin/bash

set -x

usage() {
    echo "Usage: $0 [-p <prefix> -m <manifest> -r <raft_dir> -b <bucket_basepath>]" 1>&2
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
export LENS_PROJ_ID=${projectDir}

###############################
##                           ##
## STEP THREE: RUN PIPELINE  ##
##                           ##
###############################

# run LENS
cd ${RAFT_PATH}/projects/${LENS_PROJ_ID}/logs/

PROJECT_BASEDIR_BUCKET=${b}
CONFIG_FILE_COMPLETE="${RAFT_PATH}/setup-lens-hlamajority/${projectDir}/nextflow-for-aws-lens-${LENS_PROJ_ID}-complete.config"
echo "config file complete..."
echo $CONFIG_FILE_COMPLETE

nextflow -Dnxf.pool.type=sync run "${RAFT_PATH}/projects/${LENS_PROJ_ID}/workflow/main.nf" \
                                  --project_dir ${PROJECT_BASEDIR_BUCKET} \
                                  -with-trace \
                                  -with-report \
                                  -with-dag \
                                  -with-timeline \
                                  -bucket-dir ${PROJECT_BASEDIR_BUCKET}/nextflow-workdir-${LENS_PROJ_ID} \
                                  --global_fq_dir ${PROJECT_BASEDIR_BUCKET}/fastqs \
                                  --global_bam_dir ${PROJECT_BASEDIR_BUCKET}/bams \
                                  --shared_dir ${PROJECT_BASEDIR_BUCKET}/shared \
                                  -c "${CONFIG_FILE_COMPLETE}" \
                                  --ref_dir ${PROJECT_BASEDIR_BUCKET}/references/cloud \
				  --output_dir ${PROJECT_BASEDIR_BUCKET}/projects/${LENS_PROJ_ID}/outputs \
				  --metadata_dir ${PROJECT_BASEDIR_BUCKET}/metadata \
				  -resume

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time 02-lens-nextflow-run-script.sh: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
