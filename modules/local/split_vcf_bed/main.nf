process SPLIT_VCF_BED {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools%3A1.23.1--hb2cee57_0':
        'biocontainers/bcftools:1.23.1--hb2cee57_0' }"

    input:
    tuple val(meta), path(vcf_file)
    path bed

    output:
    tuple val(meta), path ('chunk_*.vcf.gz'), emit: splited_vcfs
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
        count=0
        bcftools index -t ${vcf_file}
        cat ${bed} | awk '{FS="\t";ODS="\t"} {print \$1 ":" \$2 "-" \$3}' | while IFS= read -r line; do
            ((count+=1))
            suffix=\$(printf "%05d" \$count)
            echo "\$suffix - \$line" 1>&2
            bcftools view --regions "\$line" ${vcf_file} | bcftools norm -m- | bcftools view --drop-genotypes -Oz -o chunk_\${suffix}.vcf.gz &
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
    stub:
        """
        i="1"
        for i in \$(seq 1 10); do
            touch chunk_\${i}.vcf.gz
        done
        touch versions.yml
        """  
}