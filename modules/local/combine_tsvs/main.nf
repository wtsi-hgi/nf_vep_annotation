process COMBINE_TSVS {
    //publishDir  path: "${params.publishdir}",
    //            mode: "copy",
    //            overwrite: "true"

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/htslib%3A1.23.1--h633afcb_0':
    //    'biocontainers/htslib:1.23.1--h633afcb_0' }"

    input:
    tuple val(meta), path (vep_outputs)

    output:
    tuple val(meta), path("combined_vep_output_for_hail_qc.tsv") , emit: vep_annotations

    script:
        """
        cat ${vep_outputs.join(' ')} | sort -k1,1V -k2,2n > combined_vep_output_for_hail_qc.tsv
        """
    stub:
        """
        touch combined_vep_output_for_hail_qc.tsv
        """
}
