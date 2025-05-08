process RUN_VEP {
    //publishDir  path: "${params.publishdir}",
    //            mode: "copy",
    //            overwrite: "true"

    input:
    tuple val(meta), path (vcf_file)

    output:
    tuple val(meta), path("${vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')}"), emit: vep_vcf

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ${params.vep_cmd} ${vcf_file}
    """
}
