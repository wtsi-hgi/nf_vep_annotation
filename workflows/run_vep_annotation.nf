include { VCF_PREPROCESS } from '../modules/local/vcf_preprocess/main'
include { RUN_VEP } from '../modules/local/run_vep/main'
include { NORM_VCF; NORM_VCF as NORM_VCF_WITH_G; NO_G_VCF; SPLIT_VCF } from '../modules/local/split_vcf/main'
include {SPLIT_VCF_BED} from '../modules/local/split_vcf_bed/main'
include { BCFTOOLS_SPLIT_VEP } from '../modules/local/bcftools_split_vep/main'
include { COMBINE_TSVS } from '../modules/local/combine_tsvs/main'
include { BCFTOOLS_EXTRACT_CSQ } from '../modules/local/extract_csq/main'
include { COMBINE_CSQS } from '../modules/local/combine_csqs/main'
include { GET_CSQ_HEADER } from '../modules/local/csq_header/main'
include { BCFTOOLS_ANNOTATE } from '../modules/local/bcftools_annotation/main'




workflow RUN_VEP_ANNOTATION{
    // Create the output directory if it doesn't exist
    if (!file(params.publishdir).exists()) {
        file(params.publishdir).mkdirs()
    }
    // split VCF or not
    if("${params.split_input}"=='true'){
        vcf_file=channel.fromPath(params.vcf_infile)
        vcf_input_file=vcf_file.map{vcf_file -> [[id:'input_vcf'], vcf_file]}
        //normalise VCFs
        NO_G_VCF(vcf_input_file)
        vcf_noG=NO_G_VCF.out.no_g_vcf
        NORM_VCF(vcf_noG)
        vcf_norm=NORM_VCF.out.normolized_vcf
        //prepare normolized VCF if annotated VCF specified as an output
        if ("${params.annotate_vcf}"=='true'){
            NORM_VCF_WITH_G(vcf_file)
            norm_vcf_with_g=NORM_VCF_WITH_G.out.normolized_vcf.map{
                meta, vcf_norm -> [[id:'norm_vcf_with_genomes'], vcf_norm]
            }
        }
        //split VCF using bed or not
        if("${params.use_bed_to_split}"=='true'){
            //split VCF using bed
            bed_file=channel.fromPath(params.interval_bed)
            vcf_input = vcf_norm.map{
            vcf_norm -> [[id:'split_vcf_using_bed'], vcf_norm]
            }
            SPLIT_VCF_BED(vcf_input, bed_file)
            //prepare VCF chunks for VEP annoation
            vcf_chunks=SPLIT_VCF_BED.out.splited_vcfs.map{
            meta, splited_vcfs -> [splited_vcfs]
            }
            //ßvcf_chunks.view()
        }else{//change vcf_file here into NO_G_VCF
            //split one VCF into N chuks
            number_of_chunks=channel.value(params.number_of_chunks)
            vcf_input = number_of_chunks.combine(vcf_norm).map{
            number_of_chunks, meta, vcf_norm -> [[id:'split_vcf_to_' + number_of_chunks + '_chunks'], vcf_norm]
            }
            SPLIT_VCF(vcf_input, number_of_chunks)
            //prepare VCF chunks for VEP annoation
            vcf_chunks=SPLIT_VCF.out.splited_vcfs.map{
            meta, splited_vcfs -> [splited_vcfs]
            }
        }
        //make map for vep annotation process
        numbers=vcf_chunks.flatten().collect().map { it.size() }.map { 1..it }.flatten()
        shards=vcf_chunks.flatten().merge(numbers).map{
        vcf_chunks, numbers -> [[id:'vep_annotation_'+numbers], vcf_chunks]
        }
    }else{
        //work with VCF shards
        vcf_files = channel.fromPath("$params.vcf_in/*.vcf.gz")

        vcf_input_file = vcf_files
            .map { file ->
                def id = file.getSimpleName().replaceAll(/\.vcf$/, '')
                def meta = [
                    id: id
                ]
                tuple(meta, file)
            }

        numbers=vcf_files.collect().map { it.size() }.map { 1..it }.flatten()
        vcf_input = vcf_files.merge(numbers).map{
        vcf_file, numbers -> [[id:'vcf_preprocess_'+numbers], vcf_file]
        }
        //vcf_input_file=vcf_file.map{vcf_file -> [[id:'input_vcf'], vcf_file]}
        //normalise VCFs
        NO_G_VCF(vcf_input)
        vcf_noG=NO_G_VCF.out.no_g_vcf
        NORM_VCF(vcf_noG)
        shards=NORM_VCF.out.normolized_vcf.merge(numbers).map{
            meta, vcf_file, numbers -> [[id:'vep_annotation_'+numbers], vcf_file]
        }
        //prepare normolized VCF if annotated VCF specified as an output
        if ("${params.annotate_vcf}"=='true'){
            NORM_VCF_WITH_G(vcf_input)
            norm_vcf_with_g=NORM_VCF_WITH_G.out.normolized_vcf.merge(numbers).map{
                meta, vcf_file, numbers -> [[id:'norm_vcf_with_genomes_'+numbers], vcf_file]
            }
        }
    }
    //run VEP
    RUN_VEP(shards, params.vep_options)

    //extract VEP annotation and save as TSV files
    vep_vcfs=RUN_VEP.out.vep_vcf.merge(numbers).map{
        meta, vep_vcf, numbers -> [[id:'annotation_extraction_'+numbers], vep_vcf]
    }

    GET_CSQ_HEADER(vep_vcfs.first())
    header=GET_CSQ_HEADER.out.csq_header

    BCFTOOLS_EXTRACT_CSQ(vep_vcfs, params.transcript_mode)
    csq_tsvs=BCFTOOLS_EXTRACT_CSQ.out.vep_csq_tsv.map{
        meta, tsv_file -> [tsv_file]
    }
    reference_fasta=channel.fromPath(params.ref_fasta)
    //vep_vcfs.view()
    //vep_vcfs.combine(reference_fasta).view()
    BCFTOOLS_SPLIT_VEP(vep_vcfs.combine(reference_fasta))
    //combine VEP annotations from all shards
    tsvs=BCFTOOLS_SPLIT_VEP.out.vep_split_tsv.map{
        meta, tsf_file -> [tsf_file]
    }
    csq_tsvs=csq_tsvs.collect().map{
        tsv_files -> [[id:'annotation_concatination'], tsv_files]
    }
    vep_tsvs=tsvs.collect().map{
        tsf_files -> [[id:'annotation_concatination'], tsf_files]
    }
    COMBINE_TSVS(vep_tsvs)
    COMBINE_CSQS(csq_tsvs)
    //INDEX_TSV(COMBINE_TSVS.out.vep_annotations)
    csq=COMBINE_CSQS.out.vep_annotations
    if ("${params.annotate_vcf}"=='true'){
        //norm_vcf_with_g.view()
        //csq.view()
        //header.view()
        BCFTOOLS_ANNOTATE(norm_vcf_with_g, csq, header)
    }
    //BCFTOOLS_ANNOTATE(vcf_input_file, INDEX_TSV.out.tsv, header)
    //BCFTOOLS_ANNOTATE()
}
