process NORM_VCF {

    input:
    path vcf_file

    output:
    path "${vcf_file.baseName}.normalized.vcf.gz", emit: normolized_vcf
    path "versions.yml"                     , emit: versions

    script:
    //def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools norm --threads 10 -m- ${vcf_file} -Oz -o ${vcf_file.baseName}.normalized.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

process NO_G_VCF {

    input:
    path vcf_file

    output:
    path "${vcf_file.baseName}.noG.vcf.gz", emit: no_g_vcf
    path "versions.yml"                     , emit: versions

    script:
    //def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools view  --threads 10 --drop-genotypes ${vcf_file} -Oz -o ${vcf_file.baseName}.noG.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}


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
    def suffix_length = N.toString().length()
    """
    total_variants=\$(bcftools view --no-header ${vcf_file} | wc -l)
    chunk_size=\$(( (total_variants + ${N} - 1) / ${N} ))
    bcftools view -G -h ${vcf_file} > header.txt
    
    split_filter() { 
        { cat header.txt; cat; } | bgzip > "\$FILE"; 
    }
    export -f split_filter

    bcftools view --no-header ${vcf_file} | split --numeric-suffixes=1 --suffix-length=${suffix_length} --lines="\$chunk_size" --additional-suffix=".vcf.gz" --filter='split_filter' - chunk_

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}