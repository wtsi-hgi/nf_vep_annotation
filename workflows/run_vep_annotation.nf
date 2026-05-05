include { VCF_PREPROCESS } from '../modules/local/vcf_preprocess/main'
include { RUN_VEP } from '../modules/local/run_vep/main'
include { NORM_VCF; NO_G_VCF; SPLIT_VCF } from '../modules/local/split_vcf/main'
include {SPLIT_VCF_BED} from '../modules/local/split_vcf_bed/main'
include { BCFTOOLS_SPLIT_VEP } from '../modules/local/bcftools_split_vep/main'
include { COMBINE_TSVS } from '../modules/local/combine_tsvs/main'
include { BCFTOOLS_SPLIT_VEP_COMPLETE } from '../modules/local/bcftools_split_vep_complete/main'
include { COMBINE_TSVS_COMPLETE } from '../modules/local/combine_tsvs_complete/main'
include { PREPARE_VEP_COMPLETE_ANNOT } from '../modules/local/prepare_vep_complete_annot/main'
include { ANNOTATE_VCF_WITH_VEP_COMPLETE } from '../modules/local/annotate_vcf_with_vep_complete_csq/main'
workflow RUN_VEP_ANNOTATION{
    // Create the output directory if it doesn't exist
    if (!file(params.publishdir).exists()) {
        file(params.publishdir).mkdirs()
    }
    if (params.final_vcf_outdir && !file(params.final_vcf_outdir).exists()) {
        file(params.final_vcf_outdir).mkdirs()
    }
    // split VCF or not
    if("${params.split_input}"=='true'){
        vcf_file=channel.fromPath(params.vcf_infile)
        NORM_VCF(vcf_file)
        vcf_norm=NORM_VCF.out.normolized_vcf
        NO_G_VCF(vcf_norm)
        vcf_noG=NO_G_VCF.out.no_g_vcf
        //split VCF using bed or not
        if("${params.use_bed_to_split}"=='true'){
            //split VCF using bed
            bed_file=channel.fromPath(params.interval_bed)
            vcf_input = vcf_noG.map{
            vcf_noG -> [[id:'split_vcf_using_bed'], vcf_noG]
            }
            SPLIT_VCF_BED(vcf_input, bed_file)
            //prepare VCF chunks for VEP annoation
            vcf_chunks=SPLIT_VCF_BED.out.splited_vcfs.map{
            meta, splited_vcfs -> [splited_vcfs]
            }
            vcf_chunks.view()
        }else{//change vcf_file here into NO_G_VCF
            //split one VCF into N chuks
            number_of_chunks=channel.value(params.number_of_chunks)
            vcf_input = number_of_chunks.combine(vcf_noG).map{
            number_of_chunks, vcf_noG -> [[id:'split_vcf_to_' + number_of_chunks + '_chunks'], vcf_noG]
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
    //vep_vcfs.view()
    //vep_vcfs.combine(reference_fasta).view()
    BCFTOOLS_SPLIT_VEP(vep_vcfs.combine(reference_fasta))
    BCFTOOLS_SPLIT_VEP_COMPLETE(vep_vcfs.combine(reference_fasta))
    //combine VEP annotations from all shards
    tsvs=BCFTOOLS_SPLIT_VEP.out.vep_split_tsv.map{
        meta, tsf_file -> [tsf_file]
    }
    vep_tsvs=tsvs.collect().map{
        tsf_files -> [[id:'annotation_concatination'], tsf_files]
    }
    COMBINE_TSVS(vep_tsvs)

    //combine complete VEP annotations from all shards
    complete_tsvs=BCFTOOLS_SPLIT_VEP_COMPLETE.out.vep_complete_tsv.map{
        meta, tsf_file -> [tsf_file]
    }
    vep_complete_tsvs=complete_tsvs.collect().map{
        tsf_files -> [[id:'annotation_complete_concatination'], tsf_files]
    }
    COMBINE_TSVS_COMPLETE(vep_complete_tsvs)

    // annotate vcf_after-qc per-chromosome VCFs with VEP_COMPLETE from combined complete table
    complete_table = COMBINE_TSVS_COMPLETE.out.vep_annotations_complete.map { tsv -> [[id:'vep_complete_annot'], tsv] }
    PREPARE_VEP_COMPLETE_ANNOT(complete_table)

    vcf_after_qc_hard = Channel.fromPath("${params.final_vcf_outdir}/filtered_vcfs_combinations/*.vcf.bgz", checkIfExists: true)
        .filter { f -> !f.name.endsWith('.annotated_vep_complete.vcf.bgz') }
        .map { f -> [[id: "annotate_${f.baseName}"], "filtered_vcfs_combinations", f] }

    vcf_after_qc_stringent = Channel.fromPath("${params.final_vcf_outdir}/filtered_vcfs_combinations_stringent/*.vcf.bgz", checkIfExists: true)
        .filter { f -> !f.name.endsWith('.annotated_vep_complete.vcf.bgz') }
        .map { f -> [[id: "annotate_${f.baseName}"], "filtered_vcfs_combinations_stringent", f] }

    vcf_after_qc = vcf_after_qc_hard.mix(vcf_after_qc_stringent)

    // PREPARE emits a single tuple; turn it into a value channel so it is reused
    // for every per-chromosome VCF instead of being consumed once.
    vep_complete_annot_value = PREPARE_VEP_COMPLETE_ANNOT.out.vep_complete_annot.first()
    ANNOTATE_VCF_WITH_VEP_COMPLETE(vcf_after_qc, vep_complete_annot_value)
}
