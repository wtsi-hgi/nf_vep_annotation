include { VCF_PREPROCESS } from '../modules/local/vcf_preprocess/main'
include { RUN_VEP } from '../modules/local/run_vep/main'
include { SPLIT_VCF } from '../modules/local/split_vcf/main'
include {SPLIT_VCF_BED} from '../modules/local/split_vcf_bed/main'
include { BCFTOOLS_SPLIT_VEP } from '../modules/local/bcftools_split_vep/main'
include { COMBINE_TSVS } from '../modules/local/combine_tsvs/main'
workflow RUN_VEP_ANNOTATION{
    // Create the output directory if it doesn't exist
    if (!file(params.publishdir).exists()) {
        file(params.publishdir).mkdirs()
    }
    // split VCF or not
    if("${params.split_input}"=='true'){
        vcf_file=channel.fromPath(params.vcf_infile)
        //split VCF using bed or not
        if("${params.use_bed_to_split}"=='true'){
            //split VCF using bed
            bed_file=channel.fromPath(params.interval_bed)
            vcf_input = vcf_file.map{
            vcf_file -> [[id:'split_vcf_using_bed'], vcf_file]
            }
            SPLIT_VCF_BED(vcf_input, bed_file)
            //prepare VCF chunks for VEP annoation
            vcf_chunks=SPLIT_VCF_BED.out.splited_vcfs.map{
            meta, splited_vcfs -> [splited_vcfs]
            }
            vcf_chunks.view()
        }else{
            //split one VCF into N chuks
            number_of_chunks=channel.value(params.number_of_chunks)
            vcf_input = number_of_chunks.combine(vcf_file).map{
            number_of_chunks, vcf_file -> [[id:'split_vcf_to_' + number_of_chunks + '_chunks'], vcf_file]
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
        numbers=vcf_files.collect().map { it.size() }.map { 1..it }.flatten()
        vcf_input = vcf_files.merge(numbers).map{
        vcf_file, numbers -> [[id:'vcf_preprocess_'+numbers], vcf_file]
        }
        //normalise VCFs
        VCF_PREPROCESS(vcf_input)
        shards=VCF_PREPROCESS.out.vcf_normalised.merge(numbers).map{
        meta, vcf_file, numbers -> [[id:'vep_annotation_'+numbers], vcf_file]
        }
    }
    //run VEP
    RUN_VEP(shards, params.vep_options)

    //extract VEP annotation and save as TSV files
    vep_vcfs=RUN_VEP.out.vep_vcf.merge(numbers).map{
        meta, vep_vcf, numbers -> [[id:'annotation_extraction_'+numbers], vep_vcf]
    }
    reference_fasta=channel.fromPath(params.ref_fasta)
    BCFTOOLS_SPLIT_VEP(vep_vcfs, reference_fasta)
    //combine VEP annotations from all shards
    tsvs=BCFTOOLS_SPLIT_VEP.out.vep_split_tsv.map{
        meta, tsf_file -> [tsf_file]
    }
    vep_tsvs=tsvs.collect().map{
        tsf_files -> [[id:'annotation_concatination'], tsf_files]
    }
    COMBINE_TSVS(vep_tsvs)
}
