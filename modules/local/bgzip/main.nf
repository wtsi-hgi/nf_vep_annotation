process BGZIP {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib%3A1.23.1--h633afcb_0':
        'biocontainers/htslib:1.23.1--h633afcb_0' }"

    publishDir  path: "${params.publishdir}",
                mode: "copy",
                overwrite: "true",
                pattern: '{CSQ.tsv.gz*,WxS_QC_CSQ.tsv.gz*}'
    input:
    tuple val(meta), path (vep_tsv)

    output:
    tuple val (meta), path ("${vep_tsv}.gz"), path("${vep_tsv}.gz.tbi") , emit: vep_annotations_gziped
    script:
        """
        bgzip ${vep_tsv}
        tabix -s 1 -b 2 -e 2 ${vep_tsv}.gz
        """
    stub:
        """
        touch ${vep_tsv}.gz
        touch ${vep_tsv}.gz.tbi
        """
}
