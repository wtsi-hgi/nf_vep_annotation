#!/bin/bash

# to run: bsub < $PWD/bsub_nextflow.sh

#BSUB -o /lustre/scratch124/humgen/teams_v2/hgi/re3/DDD_WGS/farmout/%J.o
#BSUB -e /lustre/scratch124/humgen/teams_v2/hgi/re3/DDD_WGS/farmout/%J.e
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"
#BSUB -q oversubscribed
#BSUB -n 2

export HTTP_PROXY='http://wwwcache.sanger.ac.uk:3128'
export HTTPS_PROXY='http://wwwcache.sanger.ac.uk:3128'
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=22.04.0-5697

module load cellgen/nextflow/24.10.2
module load cellgen/singularity

nfdir=/path/to/nf_vep_annotation/nextflow_pipeline
workdir=/path?to/working/directory

nextflow -log nextflow.log run \
${nfdir}/main.nf \
-profile sanger \
-w $workdir \
-with-trace \
-resume
