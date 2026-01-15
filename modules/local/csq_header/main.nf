process GET_CSQ_HEADER {    
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
