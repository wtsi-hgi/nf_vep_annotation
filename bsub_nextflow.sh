#!/bin/bash

# to run: bsub < $PWD/nf_vep_annotation/bsub_nextflow.sh

#BSUB -o /nfs/users/nfs_e/eh19/work/data/gtcheck/analysis/vep_annotation_trial/logs/%J.o
#BSUB -e /nfs/users/nfs_e/eh19/work/data/gtcheck/analysis/vep_annotation_trial/logs/%J.e
#BSUB -M 8000
#BSUB -q oversubscribed
#BSUB -n 2

set -euo pipefail

REPO_ROOT="/nfs/users/nfs_e/eh19/work/data/gtcheck/analysis/vep_annotation_trial"
LOG_DIR="$REPO_ROOT/logs"
REPORT_DIR="$LOG_DIR/reports"
TRACE_DIR="$LOG_DIR/traces"
RUN_ID="${LSB_JOBID:-$(date +%s)}"

mkdir -p "$LOG_DIR" "$REPORT_DIR" "$TRACE_DIR"

export HTTP_PROXY='http://wwwcache.sanger.ac.uk:3128'
export HTTPS_PROXY='http://wwwcache.sanger.ac.uk:3128'
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=22.04.0-5697

nextflow -log "$LOG_DIR/vep-trial.$RUN_ID.log" run "$REPO_ROOT/nf_vep_annotation/main.nf" -c "$REPO_ROOT/configs/vep_trial_config.nf" -profile sanger -with-report "$REPORT_DIR/report.$RUN_ID.html" -with-trace "$TRACE_DIR/trace.$RUN_ID.txt" -resume
