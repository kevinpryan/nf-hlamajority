#!/usr/bin/env bash
set -euo pipefail
BUCKET="s3://bucket-dir/references/cloud"
DRYRUN="--dryrun"   # REMOVE after validation
#DRYRUN=""
cd ../../../../
aws s3 cp --recursive references "$BUCKET/nf-hlamajority/references" $DRYRUN
aws s3 cp assets/benchmarking_results_claeys.csv "$BUCKET/nf-hlamajority/references/majority_vote/" $DRYRUN