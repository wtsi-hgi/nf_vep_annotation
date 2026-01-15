process COMBINE_CSQS {
    publishDir  path: "${params.publishdir}",
                mode: "copy",
                overwrite: "true"

    input:
    tuple val(meta), path (vep_outputs)

    output:
    tuple val ("combined CSQ"), path ("combined_vep_csq.tsv.gz"), path("combined_vep_csq.tsv.gz.tbi") , emit: vep_annotations

    script:
        """
        cat ${vep_outputs.join(' ')} | sort -k1,1V -k2,2n > combined_vep_csq.tsv
        bgzip combined_vep_csq.tsv
        tabix combined_vep_csq.tsv.gz
        """
    stub:
        """
        touch combined_vep_csq.tsv.gz
        touch combined_vep_csq.tsv.gz.tbi
        """
}
