process COMBINE_CSQS {

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/htslib%3A1.23.1--h633afcb_0':
    //    'biocontainers/htslib:1.23.1--h633afcb_0' }"

    //publishDir  path: "${params.publishdir}",
    //            mode: "copy",
    //            overwrite: "true"
 
    input:
    tuple val(meta), path (vep_outputs)

    output:
    //tuple val ("combined CSQ"), path ("combined_vep_csq.tsv.gz"), path("combined_vep_csq.tsv.gz.tbi") , emit: vep_annotations
    tuple val(meta), path("${meta.id}.csq.tsv") , emit: vep_annotations
    script:
        """
        cat ${vep_outputs.join(' ')} | sort -k1,1V -k2,2n > ${meta.id}.csq.tsv
        """
    stub:
        """
        touch ${meta.id}.csq.tsv
        """
}
