process COMBINE_TSVS_COMPLETE {
    publishDir  path: "${params.publishdir}",
                mode: "copy",
                overwrite: "true"

    input:
    tuple val(meta), path (vep_outputs)

    output:
    path "combined_vep_output_complete.tsv", emit: vep_annotations_complete

    script:
    def first_tsv = vep_outputs[0]
    """
    head -n 1 "${first_tsv}" > combined_vep_output_complete.tsv

    for f in ${vep_outputs.join(' ')}; do
      tail -n +2 "\$f"
    done | sort -k1,1V -k2,2n >> combined_vep_output_complete.tsv
    """
}

