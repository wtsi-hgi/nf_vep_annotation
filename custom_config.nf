// Custom config for a specific run — overrides the defaults from conf/params.config.
// Used via `-c custom_config.nf` (see bsub_nextflow.sh).

params {
    input = "/path/to/input/folder/or/file"
    // Any other parameter from conf/params.config can be overridden here too, e.g.:
    // publishdir = "/path/to/results"
    // number_of_chunks = 300
    // left_align = true
    // plugins = """ \
    // --plugin SpliceRegion,Extended \
    // --plugin LoFtool,/lustre/scratch125/humgen/resources_v2/ensembl/vep/GRCh38/Plugins/LoFtool_scores.txt \
    // """
}

// This file is a full Nextflow config, not just a params file, so settings
// outside of params {} — like singularity — can be overridden here too, e.g.:
// singularity {
//     enabled    = true
//     cacheDir   = "/path/to/singularity/cache"
//     runOptions = '--bind /lustre/scratch125'
// }

// Per-process resources (cpus/memory/time/etc.) can also be overridden here
// if a run needs more than the defaults in conf/base.config, e.g. if RUN_VEP
// keeps hitting its memory limit on a larger-than-usual VCF:
// process {
//     withName: RUN_VEP {
//         memory    = { check_max( 16.GB * (task.attempt ** 1.5), 'memory' ) }
//         time      = { check_max( 96.h  * task.attempt, 'time' ) }
//         container = '/path/to/a/different/vep.sif' // override the VEP container/version if needed
//     }
// }
