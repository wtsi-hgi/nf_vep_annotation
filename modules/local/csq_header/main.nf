process GET_CSQ_HEADER {    

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools%3A1.23.1--hb2cee57_0':
        'biocontainers/bcftools:1.23.1--hb2cee57_0' }"

    input:
    tuple val(meta), path (vcf_file)
   

    output:
    tuple val("CSQ_header"), path("CSQ_header.txt"), emit: csq_header
    path "versions.yml"                     , emit: versions

    script:
        """
        bcftools view -h ${vcf_file} | grep "CSQ" > CSQ_header.txt
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
    stub:
        """
        touch CSQ_header.txt
        touch versions.yml
        """
}
