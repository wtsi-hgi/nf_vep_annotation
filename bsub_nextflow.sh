#!/bin/bash

# to run: bsub < $PWD/bsub_nextflow.sh

#BSUB -o /path/to/a/log/dir/%J.o
#BSUB -e /path/to/a/log/dir/%J.e
#BSUB -M 8000
#BSUB -q oversubscribed
#BSUB -n 2

export HTTP_PROXY='http://wwwcache.sanger.ac.uk:3128'
export HTTPS_PROXY='http://wwwcache.sanger.ac.uk:3128'
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
export NXF_VER=22.04.0-5697


nextflow run \
${pwd}/main.nf \
-profile sanger \
-with-trace \
-resume
