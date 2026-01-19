#!/bin/bash
set -x
start_time=$SECONDS
nextflow run main.nf \
	     --samplesheet samplesheet.csv \
             --references_basedir s3://path/to/bucket/with/references \
             -with-trace \
             -with-report \
             -with-dag \
             -with-timeline \
             -resume \
             -bucket-dir s3://path-to-bucket-for-workdir \
             -c /home/ec2-user/nf-hlamajority/nextflow-aws.config \
             -profile awsbatch \
	     --outdir s3://path/to/bucket/for/outdir \
	     | tee -a nextflow_run_nf-hlamajority-test.log
elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
