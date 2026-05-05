process BCFTOOLS_SPLIT_VEP_COMPLETE {
    input:
    tuple val(meta), path (vep_vcf_file), path (ref_fa)

    output:
    tuple val(meta), path("${vep_vcf_file.name.replaceAll(/\.vcf.*/, '.complete.tsv')}"), emit: vep_complete_tsv
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools norm -f ${ref_fa} -Ov ${vep_vcf_file} | \\
      bcftools +split-vep \\
        -X \\
        -d \\
        -A tab \\
        -HH \\
        -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%CSQ\\n' \\
      > ${vep_vcf_file.baseName.replaceAll(/\.vcf.*/, '.complete.tsv')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

