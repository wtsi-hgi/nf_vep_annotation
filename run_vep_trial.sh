#!/bin/bash
#BSUB -G normal
#BSUB -q "oversubscribed"
#BSUB -n 2
#BSUB -M 24G
#BSUB -R "select[mem>24G] rusage[mem=24G] span[hosts=1]"
#BSUB -o "/nfs/users/nfs_e/eh19/work/data/gtcheck/analysis/vep_annotation_trial/logs/vep-trial-%J-output.log"
#BSUB -e "/nfs/users/nfs_e/eh19/work/data/gtcheck/analysis/vep_annotation_trial/logs/vep-trial-%J-errors.log"

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PIPELINE_ROOT="$SCRIPT_DIR"
TRIAL_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
LOG_DIR="$TRIAL_ROOT/logs"
REPORT_DIR="$LOG_DIR/reports"
TRACE_DIR="$LOG_DIR/traces"
RUN_ID="${LSB_JOBID:-$(date +%s)}"

mkdir -p "$LOG_DIR" "$REPORT_DIR" "$TRACE_DIR"

echo 'load modules'
ml load cellgen/nextflow/24.10.2
ml load cellgen/singularity
ml load HGI/softpack/groups/hgi/bcftools1.21/1
ml load badger/htslib/1.22.1

echo 'specify variables'
export NXF_OPTS='-Xms6G -Xmx22G -XX:+UseSerialGC'

echo 'Starting VEP annotation trial...'
echo "Input VCF: /lustre/scratch127/humgen/projects_v2/mrc_exomes/alspac/alspacv2.joint_germline.vcf.gz"
echo "Output directory: /lustre/scratch127/humgen/teams_v2/hgi/eh19_hail_qc/vep_annotation_trial/output/"

nextflow -log "$LOG_DIR/vep-trial.${RUN_ID}.log" run "$PIPELINE_ROOT/main.nf" -c "$TRIAL_ROOT/configs/vep_trial_config.nf" -profile sanger -resume -with-report "$REPORT_DIR/report.${RUN_ID}.html" -with-trace "$TRACE_DIR/trace.${RUN_ID}.txt"
