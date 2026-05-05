process PREPARE_VEP_COMPLETE_ANNOT {
    input:
    tuple val(meta), path(complete_tsv)

    output:
    tuple val(meta), path("vep_complete_annot.tsv.gz"), path("vep_complete_annot.tsv.gz.tbi"), path("vep_complete_annot.hdr"), emit: vep_complete_annot
    path "versions.yml", emit: versions

    script:
    def outBase = (params.final_vcf_outdir ?: params.publishdir).toString().replaceAll(/\/$/, '')
    def outDir  = "${outBase}/vep_complete_annotation"
    def tmpTsv  = "${outDir}/vep_complete_annot.tsv.gz.tmp.${task.hash}"
    """
    mkdir -p "${outDir}"

    # Build a compact annotation table for bcftools annotate:
    # CHROM  POS  REF  ALT  VEP_COMPLETE
    #
    # The input complete TSV starts with:
    # CHROM POS ID REF ALT <VEP subfields...>
    #
    # We reconstruct a CSQ-like string by joining all subfields (col 6..NF) with '|'.
    awk -F'\\t' 'BEGIN{OFS=\"\\t\"} NR==1{next} { csq=\$6; for(i=7;i<=NF;i++) csq=csq \"|\" \$i; print \$1,\$2,\$4,\$5,csq }' ${complete_tsv} \\
      | bgzip -c > "${tmpTsv}"

    mv -f "${tmpTsv}" "${outDir}/vep_complete_annot.tsv.gz"
    tabix -s1 -b2 -e2 "${outDir}/vep_complete_annot.tsv.gz"

    cat > "${outDir}/vep_complete_annot.hdr" <<'EOF'
##INFO=<ID=VEP_COMPLETE,Number=.,Type=String,Description="VEP annotations (all consequences) reconstructed from combined_vep_output_complete.tsv; per transcript records separated by commas, subfields by |">
EOF

    # Symlink into work dir so Nextflow can track outputs without duplicating data
    ln -sf "${outDir}/vep_complete_annot.tsv.gz" vep_complete_annot.tsv.gz
    ln -sf "${outDir}/vep_complete_annot.tsv.gz.tbi" vep_complete_annot.tsv.gz.tbi
    ln -sf "${outDir}/vep_complete_annot.hdr" vep_complete_annot.hdr

    cat > versions.yml <<'END_VERSIONS'
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //' | cut -d' ' -f1)
END_VERSIONS
    """
}

