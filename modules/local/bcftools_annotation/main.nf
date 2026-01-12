process BCFTOOLS_ANNOTATE {    
    input:
    tuple val(meta), path (vcf_file)
    tuple val(meta), path (vep_tsv), path (vep_index)
    tuple path(header)
   

    output:
    tuple val(meta), path("${vep_vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')}"), path(vcf_index), emit: annotated_vcf
    path "versions.yml"                     , emit: versions

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def vcf_index = ${vep_vcf_file.name.replaceAll(/\.vcf.gz/, '.vep.vcf.gz.tbi')}
        """
        bcftools norm -f ${ref_fa} -Ov ${vep_vcf_file} | bcftools +split-vep -s worst -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CSQ\t%Consequence\t%SYMBOL\t%HGNC_ID\n' > ${vep_vcf_file.baseName.replaceAll(/\.vcf.*/, '.tsv')}
        bcftools annotate -a ${vcf_file} -h ${header} -c CHROM,POS,REF,ALT,INFO/CSQ -Oz -o ${vep_vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')} ${vcf_file}
        bcftools index -t ${vcf_file}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
    stub:
        """
        touch ${vep_vcf_file.name.replaceAll(/\.vcf/, '.vep.vcf')}
        touch ${vcf_index}
        touch versions.yml
        """
}
