process VCF_PREPROCESS {

    input:
    tuple val(meta), path (vcf_file)

    output:
    tuple val(meta), path("${vcf_file.name.replaceAll(/\.vcf/, '.normalised.vcf')}"), emit: vcf_normalised
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools norm -m- ${vcf_file} | bcftools view -G -Oz -o ${vcf_file.name.replaceAll(/\.vcf/, '.normalised.vcf')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
