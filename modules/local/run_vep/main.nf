process RUN_VEP {

    input:
    tuple val(meta), path (vcf_file)
    val vep_options

    output:
    tuple val(meta), path("${vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')}"), emit: vep_vcf

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def vcf_out="${vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')}"
    """
    vep \
    --cache \
    ${vep_options} \
    --offline \
    --format vcf \
    -i ${vcf_file} \
    --fork 4 \
    --everything \
    --vcf \
    -o ${vcf_out} \
    --compress_output bgzip \
    --allele_number \
    --verbose
    """
}
