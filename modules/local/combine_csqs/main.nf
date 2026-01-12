process COMBINE_CSQS {
    publishDir  path: "${params.publishdir}",
                mode: "copy",
                overwrite: "true"

    input:
    tuple val(meta), path (vep_outputs)

    output:
    path "combined_vep_csq.tsv" , emit: vep_annotations

    script:
        """
        cat ${vep_outputs.join(' ')} | sort -k1,1V -k2,2n > combined_vep_csq.tsv
        """
    stub:
        """
        touch combined_vep_csq.tsv
        """
}
