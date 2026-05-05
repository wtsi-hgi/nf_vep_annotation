nextflow.enable.dsl=2

//include { RUN_GLIMPSE } from './workflows/run_glimpse'
include { RUN_VEP_ANNOTATION } from './workflows/run_vep_annotation'
include { PREPARE_VEP_COMPLETE_ANNOT } from './modules/local/prepare_vep_complete_annot/main'
include { ANNOTATE_VCF_WITH_VEP_COMPLETE } from './modules/local/annotate_vcf_with_vep_complete_csq/main'

workflow MAIN {
    //RUN_GLIMPSE ()
    RUN_VEP_ANNOTATION ()
}

/*
 * Run only the PREPARE_VEP_COMPLETE_ANNOT step.
 *
 * Usage:
 *   nextflow run ./nf_vep_annotation/main.nf -c configs/vep_trial_config.nf -profile sanger -entry PREPARE_ONLY -resume
 */
workflow PREPARE_ONLY {
    complete_tsv = Channel.fromPath("${params.publishdir}/combined_vep_output_complete.tsv", checkIfExists: true)
        .map { f -> [[id:'vep_complete_annot'], f] }
    PREPARE_VEP_COMPLETE_ANNOT(complete_tsv)
}

/*
 * Prepare the annotation table and annotate all per-chromosome VCFs in
 * ${params.final_vcf_outdir}/filtered_vcfs_combinations*.
 *
 * Outputs are written back into the same folders with suffix:
 *   *.annotated_vep_complete.vcf.bgz (+ .tbi)
 *
 * Usage:
 *   nextflow run ./nf_vep_annotation/main.nf -c configs/vep_trial_config.nf -profile sanger -entry ANNOTATE_ONLY -resume
 */
workflow ANNOTATE_ONLY {
    complete_tsv = Channel.fromPath("${params.publishdir}/combined_vep_output_complete.tsv", checkIfExists: true)
        .map { f -> [[id:'vep_complete_annot'], f] }
    PREPARE_VEP_COMPLETE_ANNOT(complete_tsv)

    vcf_after_qc_hard = Channel.fromPath("${params.final_vcf_outdir}/filtered_vcfs_combinations/*.vcf.bgz", checkIfExists: true)
        .filter { f -> !f.name.endsWith('.annotated_vep_complete.vcf.bgz') }
        .map { f -> [[id: "annotate_${f.baseName}"], "filtered_vcfs_combinations", f] }

    vcf_after_qc_stringent = Channel.fromPath("${params.final_vcf_outdir}/filtered_vcfs_combinations_stringent/*.vcf.bgz", checkIfExists: true)
        .filter { f -> !f.name.endsWith('.annotated_vep_complete.vcf.bgz') }
        .map { f -> [[id: "annotate_${f.baseName}"], "filtered_vcfs_combinations_stringent", f] }

    vcf_after_qc = vcf_after_qc_hard.mix(vcf_after_qc_stringent)

    // PREPARE emits one tuple; make it a value channel so all VCFs are annotated.
    vep_complete_annot_value = PREPARE_VEP_COMPLETE_ANNOT.out.vep_complete_annot.first()
    ANNOTATE_VCF_WITH_VEP_COMPLETE(vcf_after_qc, vep_complete_annot_value)
}

workflow {
    MAIN ()
}

