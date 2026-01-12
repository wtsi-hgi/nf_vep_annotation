process BCFTOOLS_EXTRACT_CSQ {    
    input:
    tuple val(meta), path (vep_vcf_file)
   

    output:
    tuple val(meta), path("${vep_vcf_file.name.replaceAll(/\.vcf.*/, '.csq.tsv')}"), emit: vep_csq_tsv
    path "versions.yml"                     , emit: versions

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

        """    
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\n' ${vep_vcf_file} > ${vep_vcf_file.baseName.replaceAll(/\.vcf.*/, '.csq.tsv')}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
    stub:
        """
        touch ${vep_vcf_file.name.replaceAll(/\.vcf.*/, '.csq.tsv')}
        touch versions.yml
        """
}
