nextflow.enable.dsl=2

//include { RUN_GLIMPSE } from './workflows/run_glimpse'
include { RUN_VEP_ANNOTATION } from './workflows/run_vep_annotation'

workflow MAIN {
    //RUN_GLIMPSE ()
    RUN_VEP_ANNOTATION ()
}

workflow {
    MAIN ()
}

