#!/bin/bash

# to run: bsub < $PWD/bsub_nextflow.sh

#BSUB -o "VEP-%J-output.log"
#BSUB -e "VEP-%J-errors.log"
#BSUB -q "long"
#BSUB -n 2
#BSUB -M 8G
#BSUB -R "select[mem>8G] rusage[mem=8G]"

ml load cellgen/nextflow/24.10.2
ml load cellgen/singularity

export NXF_OPTS='-Xms6G -Xmx22G -XX:+UseSerialGC'

nfdir=/path/to/nf_vep_annotation/nextflow_pipeline
workdir=/path/to/working/directory

nextflow \
	run ${nfdir}/main.nf \
	-profile sanger \
	-c custom_config.nf \
	-w $workdir \
	-resume \
	-log vep.${LSB_JOBINDEX}.log \
	-with-report vep_report.${LSB_JOBID}.html \
	-with-trace \
