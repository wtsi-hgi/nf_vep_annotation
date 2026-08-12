process NORM_VCF {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools%3A1.23.1--hb2cee57_0':
        'biocontainers/bcftools:1.23.1--hb2cee57_0' }"

    input:
    tuple val(meta), path (vcf_file)

    output:
    tuple val(meta), path ("${vcf_file.baseName}.normalized.vcf.gz"), emit: normolized_vcf
    path "versions.yml"                     , emit: versions

    script:
        //def args = task.ext.args ?: ''
        //def prefix = task.ext.prefix ?: "${meta.id}"

        """
        bcftools norm --threads ${task.cpus} -m- ${vcf_file} -Oz -o ${vcf_file.baseName}.normalized.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
    stub:
        """
        touch ${vcf_file.baseName}.normalized.vcf.gz
        touch versions.yml
        """  
}


process NORM_VCF_LEFT_ALIGN {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools%3A1.23.1--hb2cee57_0':
        'biocontainers/bcftools:1.23.1--hb2cee57_0' }"

    input:
    tuple val(meta), path (vcf_file), path (reference_fasta)

    output:
    tuple val(meta), path("${vcf_file.baseName}.normalized.vcf.gz"), emit: la_vcf
    path "versions.yml", emit: versions

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

        """    
        bcftools norm --threads ${task.cpus} -m- -f ${reference_fasta} ${vcf_file} -Oz -o ${vcf_file.baseName}.normalized.vcf.gz
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
    stub:
        """
        touch ${vcf_file.baseName}.normalized.vcf.gz
        touch versions.yml
        """
}


process NO_G_VCF {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools%3A1.23.1--hb2cee57_0':
        'biocontainers/bcftools:1.23.1--hb2cee57_0' }"

    input:
    tuple val(meta), path(vcf_file)

    output:
    tuple val(meta), path ("${vcf_file.baseName}.noG.vcf.gz"), emit: no_g_vcf
    path "versions.yml"                     , emit: versions

    script:
        //def args = task.ext.args ?: ''
        //def prefix = task.ext.prefix ?: "${meta.id}"

        """
        bcftools view  --threads ${task.cpus} --drop-genotypes ${vcf_file} -Oz -o ${vcf_file.baseName}.noG.vcf.gz
    
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
    stub:
        """
        touch ${vcf_file.baseName}.noG.vcf.gz
        touch versions.yml
        """  
}


process SPLIT_VCF {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools%3A1.23.1--hb2cee57_0':
        'biocontainers/bcftools:1.23.1--hb2cee57_0' }"

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

        #bcftools view --no-header ${vcf_file} | split --numeric-suffixes=1 --suffix-length=${suffix_length} --lines="\$chunk_size" --additional-suffix=".vcf.gz" --filter='split_filter' - chunk_

        bcftools view --no-header ${vcf_file} | \
        awk -v chunk_size="\$chunk_size" '
        BEGIN {
            chunk = 1
            line_count = 0
        }
        {
            if (line_count % chunk_size == 0) {
                if (cmd)
                    close(cmd)

                filename = sprintf("chunk_%0'${suffix_length}'d.vcf.gz", chunk)
                cmd = "cat header.txt - | bgzip > " filename
                chunk++
            }

            print | cmd
            line_count++
        }
        END {
            if (cmd)
                close(cmd)
        }'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
    stub:
        """
        i="1"
        for i in \$(seq 1 ${N}); do
            touch chunk_\${i}.vcf.gz
        done
        touch versions.yml
        """  
}