include { VCF_PREPROCESS } from '../modules/local/vcf_preprocess/main'
include { RUN_VEP } from '../modules/local/run_vep/main'
include { NORM_VCF; NORM_VCF as NORM_VCF_WITH_G; NORM_VCF_LEFT_ALIGN; NORM_VCF_LEFT_ALIGN as NORM_VCF_WITH_G_LEFT_ALIGN; NO_G_VCF; SPLIT_VCF } from '../modules/local/split_vcf/main'
include {SPLIT_VCF_BED} from '../modules/local/split_vcf_bed/main'
include { BCFTOOLS_SPLIT_VEP } from '../modules/local/bcftools_split_vep/main'
include { COMBINE_TSVS } from '../modules/local/combine_tsvs/main'
include { BCFTOOLS_EXTRACT_CSQ } from '../modules/local/extract_csq/main'
include { COMBINE_CSQS; COMBINE_CSQS as COMBINE_ALL_CSQS } from '../modules/local/combine_csqs/main'
include { GET_CSQ_HEADER } from '../modules/local/csq_header/main'
include { BCFTOOLS_ANNOTATE } from '../modules/local/bcftools_annotation/main'
include {BGZIP; BGZIP as BGZIP2; BGZIP as BGZIP3} from '../modules/local/bgzip/main'
workflow RUN_VEP_ANNOTATION{
    // Create the output directory if it doesn't exist
    if (!file(params.publishdir).exists()) {
        file(params.publishdir).mkdirs()
    }

    if (params.hail_tsv && !params.left_align) {
        log.error "To generate Hail TSV, left alignment is required!"
    }
    if (params.split_input && params.number_of_chunks < 2) {
        log.error "To split VCF into chunks, number_of_chunks must be greater than 1!"
    }
    if (params.left_align && (!file(params.ref_fasta).exists() || !file(params.ref_fasta).isFile())) {
        log.error "Reference FASTA file is required for left alignment!"
    }

    if (file(params.input).exists()){
        if (file(params.input).isFile()){//one VCF file as an input
            vcf_input=channel.fromPath(params.input).map{vcf_file -> [[id:'input_vcf'], vcf_file]}
        }else if (file(params.input).isDirectory()){//multiple VCF files as an input
            vcf_input = channel.fromPath("${params.input}/*.{vcf,vcf.gz,vcf.bgz,bcf,bcf.gz,bcf.bgz}")
                .map { vcf_file ->
                    def id = vcf_file.getName().replaceAll(/\.(vcf|bcf)(\.gz|\.bgz)?$/, '')
                    [[id: id], vcf_file]
            }
        }else{
            log.error "Input path is neither a file nor a directory: ${params.input}"
        }
        //remove genotypes from VCF files and normalize VCF files
        NO_G_VCF(vcf_input)
        vcf_noG=NO_G_VCF.out.no_g_vcf

        //normalize vcf files with or without left align indels
        reference_fasta=channel.fromPath(params.ref_fasta)
        if (params.left_align) {//left align indels in VEP annotation. Required for hail file
            NORM_VCF_LEFT_ALIGN(vcf_noG.combine(reference_fasta))
            vcf_norm=NORM_VCF_LEFT_ALIGN.out.la_vcf
        } else {//normalize VCF files without left align indels
            NORM_VCF(vcf_noG)
            vcf_norm=NORM_VCF.out.normolized_vcf
        }

        //prepare normolized VCF if annotated VCF specified as an output
        if (params.annotate_vcf){
            if (params.left_align) {//left align indels in VEP annotated VCF files
                reference_fasta=channel.fromPath(params.ref_fasta)
                NORM_VCF_WITH_G_LEFT_ALIGN(vcf_input.combine(reference_fasta))
                norm_vcf_with_g=NORM_VCF_WITH_G_LEFT_ALIGN.out.la_vcf
            }else{//normalize VCF files without left align indels
                NORM_VCF_WITH_G(vcf_input)
                norm_vcf_with_g=NORM_VCF_WITH_G.out.normolized_vcf
            }
        }

        //split VCF files into N chunks for faster VEP annotation
        if(params.split_input){
            number_of_chunks=channel.value(params.number_of_chunks)
            SPLIT_VCF(vcf_norm, number_of_chunks)
            vcf_chunks=SPLIT_VCF.out.splited_vcfs.flatMap { meta, files ->
                files.collect { file -> [meta, file] }
            }
            vcf_channel=vcf_chunks
        }else{
            vcf_channel=vcf_norm
        }
        //run VEP
        def vep_options = """--dir_cache ${params.vep_data_dir} \
            --assembly ${params.assembly} \
            --fasta ${params.vep_fasta} \
            --dir_plugins ${params.vep_plugins_dir} ${params.plugins_to_use}"""
        RUN_VEP(vcf_channel, vep_options)
        vep_vcfs=RUN_VEP.out.vep_vcf

        //extract CSQ header from the first VEP annotated VCF file
        if (params.annotate_vcf || params.hail_tsv){
            GET_CSQ_HEADER(vep_vcfs.first())
            header=GET_CSQ_HEADER.out.csq_header
        }

        //extract VEP annotation and save as a TSV file
        BCFTOOLS_EXTRACT_CSQ(vep_vcfs, params.transcript_mode)
        //combine by meta id. turned off
        csq_tsvs=BCFTOOLS_EXTRACT_CSQ.out.vep_csq_tsv.groupTuple()

        //make VEP annotation as a new INFO field in the original VCF file
        if (params.annotate_vcf){
            if(params.split_input){
                COMBINE_CSQS(csq_tsvs)
                BGZIP(COMBINE_CSQS.out.vep_annotations)
            }else{
                BGZIP(csq_tsvs)
            }
            csq=BGZIP.out.vep_annotations_gziped
            norm_vcf_with_g_and_csq=norm_vcf_with_g.combine(csq, by: 0)
            BCFTOOLS_ANNOTATE(norm_vcf_with_g_and_csq, header)
        }

        if (params.csq_tsv){
            all_csq_tsvs=BCFTOOLS_EXTRACT_CSQ.out.vep_csq_tsv.map{
                meta, tsv_file -> [tsv_file]
            }
            all_csq_tsvs=all_csq_tsvs.collect().map{
                tsv_files -> [[id:'annotation_concatination'], tsv_files]
            }
            COMBINE_ALL_CSQS(all_csq_tsvs)
            all_csqs = COMBINE_ALL_CSQS.out.vep_annotations.map { meta, file ->
                def target = file.resolveSibling('CSQ.tsv')
                file.copyTo(target)
                [meta, target]
            }
            BGZIP2(all_csqs)
        }

        //make tsv files for hail qc
        if (params.hail_tsv && params.left_align){
            BCFTOOLS_SPLIT_VEP(vep_vcfs)
            //combine VEP annotations from all shards
            tsvs=BCFTOOLS_SPLIT_VEP.out.vep_split_tsv.map{
                meta, tsf_file -> [tsf_file]
            }
            vep_tsvs=tsvs.collect().map{
                tsf_files -> [[id:'annotation_concatination'], tsf_files]
            }
            if (file(params.input).isDirectory() || params.split_input){
                COMBINE_TSVS(vep_tsvs)
                hail_csq_tsv=COMBINE_TSVS.out.vep_annotations
            }else{
                hail_csq_tsv=vep_tsvs.map { meta, file ->
                    def target = file.resolveSibling('WxS_QC_CSQ.tsv')
                    file.copyTo(target)
                    [meta, target]
                }
            }
            BGZIP3(hail_csq_tsv)
        }
    }else{
        log.error "Input path doesn't exist: ${params.input}"
    }
}
