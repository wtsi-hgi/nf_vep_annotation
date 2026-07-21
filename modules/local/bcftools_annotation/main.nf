process BCFTOOLS_ANNOTATE {
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools%3A1.23.1--hb2cee57_0':
        'biocontainers/bcftools:1.23.1--hb2cee57_0' }"

    publishDir  path: "${params.publishdir}",
                mode: "copy",
                pattern: '*.vep.vcf.gz, *.vep.vcf.gz.tbi',
                overwrite: "true"

    input:
    tuple val(meta), path (vcf_file)
    tuple val(meta), path (vep_tsv), path (vep_index)
    tuple val(meta), path(header)
   

    output:
    tuple val(meta), path("${vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')}"), path("${vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')}.tbi"), emit: annotated_vcf
    path "versions.yml"                     , emit: versions

    script:
        def args = task.ext.args ?: ''
        //def prefix = task.ext.prefix ?: "${meta.id}"
        //def vcf_index = ${vcf_file.name.replaceAll(/\.vcf.gz/, '.vep.vcf.gz.tbi')}
        """
        bcftools annotate -a ${vep_tsv} -h ${header} -c CHROM,POS,REF,ALT,INFO/CSQ -Oz -o ${vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')} ${vcf_file}
        bcftools index -t ${vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
    stub:
        """
        touch ${vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')}
        touch ${vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')}.tbi
        touch versions.yml
        """
}
