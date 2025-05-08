process BCFTOOLS_ALL_VEP {
    input:
    tuple val(meta), path (vep_vcf_file)
    path ref_fa

    output:
    tuple val(meta), path("${vep_vcf_file.name.replaceAll(/\.vcf.*/, '.tsv')}"), emit: vep_all_tsv
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools norm -f ${ref_fa} -Ov ${vep_vcf_file} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\n' > ${vep_vcf_file.baseName.replaceAll(/\.vcf.*/, '.tsv')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
