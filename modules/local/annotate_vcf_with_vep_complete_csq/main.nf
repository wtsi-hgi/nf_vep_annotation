process ANNOTATE_VCF_WITH_VEP_COMPLETE {
    input:
    tuple val(meta), val(out_subdir), path(vcf_in)
    tuple val(meta2), path(annot_tsv_gz), path(annot_tsv_tbi), path(annot_hdr)

    output:
    tuple val(meta), path("${out_subdir}/${vcf_in.name.replaceAll(/\.vcf\.(bgz|gz)$/, '.annotated_vep_complete.vcf.bgz')}"), emit: annotated_vcf
    path("${out_subdir}/${vcf_in.name.replaceAll(/\.vcf\.(bgz|gz)$/, '.annotated_vep_complete.vcf.bgz.tbi')}"), emit: annotated_vcf_index
    path "versions.yml", emit: versions

    script:
    def out_vcf = vcf_in.name.replaceAll(/\.vcf\.(bgz|gz)$/, '.annotated_vep_complete.vcf.bgz')
    def outBase = (params.final_vcf_outdir ?: params.publishdir).toString().replaceAll(/\/$/, '')
    def outDir  = "${outBase}/${out_subdir}"
    def outPath = "${outDir}/${out_vcf}"
    def tmpPath = "${outDir}/${out_vcf}.tmp.${task.hash}"
    """
    mkdir -p "${outDir}"
    mkdir -p "${out_subdir}"

    # Overwrite INFO/CSQ using the prepared annotation table.
    # Keep output as a separate file (*.annotated_vep_complete.vcf.bgz).
    # NOTE: We intentionally do not add a new CSQ header here (-h) because
    # input VCFs already contain the Ensembl-style CSQ header definition.
    bcftools annotate \\
      -a ${annot_tsv_gz} \\
      -c CHROM,POS,REF,ALT,INFO/CSQ \\
      -Oz \\
      -o "${tmpPath}" \\
      ${vcf_in}

    mv -f "${tmpPath}" "${outPath}"
    tabix -p vcf "${outPath}"

    # Symlink into work dir so Nextflow can track outputs without duplicating data
    ln -sf "${outPath}" "${out_subdir}/${out_vcf}"
    ln -sf "${outPath}.tbi" "${out_subdir}/${out_vcf}.tbi"

    cat > versions.yml <<'END_VERSIONS'
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //' | cut -d' ' -f1)
END_VERSIONS
    """
}

