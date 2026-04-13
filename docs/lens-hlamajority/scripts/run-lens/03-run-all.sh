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

echo "Running 01-lens-setup-script.sh..."
bash 01-lens-setup-script.sh -p "${p}" \
	              -m "${m}" \
		      -r "${r}" \
		      -b "${b}" \
		      || { echo "ERROR: 01-lens-setup-script.sh failed with exit code $?. Aborting."; exit 1; }

echo "01-lens-setup-script.sh complete. Running 02-lens-nextflow-run-script.sh..."

bash 02-lens-nextflow-run-script.sh -p "${p}" \
                      -m "${m}" \
                      -r "${r}" \
                      -b "${b}" \
		      || { echo "ERROR: 02-lens-nextflow-run-script.sh failed with exit code $?. Aborting."; exit 1; }

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time 03-run-all.sh: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

