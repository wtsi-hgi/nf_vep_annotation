process SPLIT_VCF {

    input:
    tuple val(meta), path(vcf_file)
    val N

    output:
    tuple val(meta), path ('chunk_*.vcf.gz'), emit: splited_vcfs
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    total_variants=\$(bcftools norm -m- ${vcf_file} | bcftools view --no-header --drop-genotypes | wc -l)
    chunk_size=\$(( (total_variants + ${N} - 1) / ${N} ))
    bcftools view -G -h ${vcf_file} > header.txt
    
    split_filter() { 
        { cat header.txt; cat; } | bgzip > "\$FILE"; 
    }
    export -f split_filter

    bcftools norm -m- ${vcf_file} | bcftools view --no-header --drop-genotypes | split --numeric-suffixes=1 --lines="\$chunk_size" --additional-suffix=".vcf.gz" --filter='split_filter' - chunk_

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
