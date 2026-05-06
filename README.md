# nf_vep_annotation
Run VEP and extract variant annotations.

## Pipeline overview
The pipeline takes either one large input VCF or a directory of pre-sharded VCFs, runs Ensembl VEP on each shard, extracts tabular annotations, and combines them into dataset-level outputs. In this branch it also prepares a "complete" VEP annotation table and writes annotated per-chromosome VCFs back to the final output directories.

## Steps
1. Input handling:
   Accept either a single VCF (`params.vcf_infile`) or a directory of shard VCFs (`params.vcf_in`). If `split_input` is true, the input VCF is prepared for splitting; otherwise existing shards are processed directly.
2. VCF normalization and genotype dropping:
   `NORM_VCF` normalizes multiallelic records with `bcftools norm`. `NO_G_VCF` removes genotype columns with `bcftools view --drop-genotypes` so the annotation path is lighter.
3. Sharding:
   `SPLIT_VCF` splits a normalized VCF into a chosen number of chunks, or `SPLIT_VCF_BED` splits by BED intervals when `use_bed_to_split` is enabled.
4. VEP annotation:
   `RUN_VEP` runs Ensembl VEP in offline mode on each shard and writes compressed annotated VCFs.
5. Worst-consequence table extraction:
   `BCFTOOLS_SPLIT_VEP` uses `bcftools +split-vep -s worst` to convert each VEP-annotated VCF into a TSV with one selected consequence per variant/transcript record.
6. Combined summary table:
   `COMBINE_TSVS` concatenates and coordinate-sorts the shard TSVs into `combined_vep_output.tsv` in `params.publishdir`.
7. Complete annotation table:
   `BCFTOOLS_SPLIT_VEP_COMPLETE` and `COMBINE_TSVS_COMPLETE` build a second TSV that keeps the full VEP consequence content rather than only the reduced summary view.
8. Prepare annotation files for VCF rewriting:
   `PREPARE_VEP_COMPLETE_ANNOT` converts the combined complete TSV into a bgzipped/tabixed annotation table plus header metadata for `bcftools annotate`.
9. Rewrite final per-chromosome VCFs:
   `ANNOTATE_VCF_WITH_VEP_COMPLETE` applies the prepared table to VCFs under `params.final_vcf_outdir`, producing `*.annotated_vep_complete.vcf.bgz` files and indexes.

## Entry workflows
- `MAIN`: runs the full annotation pipeline (Steps 1-9).
- `PREPARE_ONLY`: runs Step 8 only, rebuilding the prepared complete-annotation table from `combined_vep_output_complete.tsv`.
- `ANNOTATE_ONLY`: runs Step 9 only, applying the prepared complete annotation to the final per-chromosome VCF outputs.

## Trial entrypoint script
- `run_vep_trial.sh` is the cluster launcher for the trial run. It is an LSF batch script (`#BSUB` directives) that requests 2 cores and 24 GB RAM.
- The script loads the required Nextflow, Singularity, `bcftools`, and `htslib` modules, creates `logs/`, `logs/reports/`, and `logs/traces/` directories under the trial root, and sets `NXF_OPTS` for the Nextflow JVM.
- It runs `nextflow run main.nf` with `-c "$TRIAL_ROOT/configs/vep_trial_config.nf"`, `-profile sanger`, and `-resume`, and also writes a Nextflow log, HTML report, and trace file for the run.
- Because it launches `main.nf` without specifying another entry workflow, this script uses the `MAIN` entry workflow by default (Steps 1-9).

Note: the `dataset-alspac` branch is kept as a dataset-specific record of ALSPAC-oriented annotation changes that generated the results in `/lustre/scratch127/humgen/projects_v2/mrc_exomes/alspac/vep_new_run`.
