# nf_vep_annotation

A Nextflow (DSL2) pipeline that annotates VCF files with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), with optional support for splitting large VCFs for parallel annotation, re-injecting VEP's `CSQ` annotation back into the original VCF, and producing flattened TSV outputs (including a TSV file for [WxS-QC](https://github.com/wtsi-hgi/wxs-qc)).

## Requirements

- [Nextflow](https://www.nextflow.io/) (tested with `24.10.2`)
- [Singularity](https://sylabs.io/singularity/) (all processes run in containers; Docker also works if you adjust `nextflow.config`)
- Access to the reference/resource paths used by VEP (cache, plugins, FASTA) and, if left-aligning, a reference FASTA for `bcftools norm`

## Usage

### Quick start

```bash
nextflow run main.nf -profile sanger --input /path/to/your.vcf.gz
```

Or point `--input` at a folder of VCFs to annotate them all in one run:

```bash
nextflow run main.nf -profile sanger --input /path/to/vcf_folder/
```

For a real run, you'll typically want your own config file layered on top of the defaults rather than passing every flag on the command line.

### Customizing a run

Rather than editing `conf/params.config` directly, create your own config file and pass it with `-c`:

```bash
nextflow run main.nf -profile sanger -c custom_config.nf
```

A template is provided at `custom_config.nf`. Since it's a full Nextflow config (not just a params file), you can override:

- **Any parameter**, inside `params { }` — including which VEP plugins to run:
  ```groovy
  params {
      input            = "/path/to/your/input"
      left_align       = true
      number_of_chunks = 300
      plugins          = "--plugin LoFtool,/lustre/scratch125/humgen/resources_v2/ensembl/vep/GRCh38/Plugins/LoFtool_scores.txt"
  }
  ```
  See [VEP plugins](#vep-plugins) below for more ready-to-use plugin examples you can paste in here.
- **Singularity settings**, e.g. the container cache location:
  ```groovy
  singularity {
      cacheDir = "/path/to/singularity/cache"
  }
  ```
- **Per-process resources or containers**, if a run needs more than the defaults in `conf/base.config` — e.g. if `RUN_VEP` runs out of memory on a larger-than-usual VCF:
  ```groovy
  process {
      withName: RUN_VEP {
          memory    = { check_max( 16.GB * (task.attempt ** 1.5), 'memory' ) }
          container = '/path/to/a/different/vep.sif'
      }
  }
  ```

### Running on LSF

`bsub_nextflow.sh` is an example submission script for an LSF cluster (`bsub < bsub_nextflow.sh`), which loads Nextflow/Singularity modules and runs with `-profile sanger -c custom_config.nf -resume`. Adjust the paths and `#BSUB` resource requests at the top of the script for your own job.

## Parameters

### Input

`--input` accepts either:

- **A single VCF/BCF file** — annotated on its own.
- **A folder** — every file in it matching `.vcf`, `.vcf.gz`, `.vcf.bgz`, `.bcf`, `.bcf.gz`, or `.bcf.bgz` is picked up and annotated independently. Anything else in that folder is ignored.

There is no built-in default for `--input` — you must supply one, or the pipeline will fail with "Input path doesn't exist".

### Run modes

These are independent toggles in `conf/params.config` — enable any combination you need.

| Flag | Default | Effect |
|---|---|---|
| `annotate_vcf` | `true` | Produce a VEP-annotated VCF: the original VCF (with genotypes) plus a `CSQ` INFO field, bgzipped and indexed. |
| `csq_tsv` | `false` | Produce a single combined TSV of every variant's CSQ annotation, across all input files/chunks. |
| `wxs_tsv` | `false` | Produce a TSV formatted for the WxS-QC pipeline. **Requires `left_align = true`.** |
| `left_align` | `false` | Left-align indels (via `bcftools norm -f ref_fasta`) before annotating. Required for `wxs_tsv`. Requires `ref_fasta` to point to a real file. |
| `split_input` | `false` | Split each input VCF into `number_of_chunks` pieces before annotating, for parallelism. |

#### Outputs

Published to `publishdir` (default `${launchDir}/results`):

| Mode | Files published |
|---|---|
| `annotate_vcf = true` | `<name>.vep.vcf.gz` + `.tbi` — one per input file |
| `csq_tsv = true` | `CSQ.tsv.gz` + `.tbi` — one combined file across all inputs |
| `wxs_tsv = true` (+ `left_align = true`) | `WxS_QC_CSQ.tsv.gz` + `.tbi` — one combined file across all inputs |

### Parameters reference

All defaults live in `conf/params.config`. The ones you'll most commonly want to override:

| Parameter | Purpose |
|---|---|
| `input` | VCF file or folder (see [Input](#input)) |
| `publishdir` | Where final results are copied to |
| `split_input` / `number_of_chunks` | Parallelise annotation of a large VCF |
| `transcript_mode` | `all` \| `worst` \| `primary` \| `mane` — passed as the `-s` (`--select`) option to `bcftools +split-vep` when extracting `CSQ` into TSV, controlling which transcript(s) are reported per variant |
| `annotate_vcf` / `csq_tsv` / `wxs_tsv` / `left_align` | Run modes, see above |
| `ref_fasta` | Reference FASTA for left-alignment (only read when `left_align = true`) |
| `vep_data_dir` / `vep_plugins_dir` / `vep_fasta` / `assembly` | VEP cache/plugin/reference locations |
| `plugins` | VEP plugin/custom-annotation flags actually passed to VEP (empty by default — see below) |

Note: `workdir` in `conf/params.config` is **informational only** — it does not set Nextflow's actual work directory. Use the `-w <path>` CLI flag for that.

#### VEP plugins

`plugins` is the value actually passed to VEP's `--dir_plugins` invocation, and is empty by default (no plugins run). `plugins_examples` in the same file is a **reference only** — the pipeline never reads it. It's a ready-made menu of correctly-formatted `--plugin ...` / `--custom ...` invocations (CADD, LoF/LOFTEE, REVEL, SpliceAI, dbNSFP, ClinVar, and more) for this project's resource paths. To enable a plugin, copy the relevant line(s) out of `plugins_examples` and paste them into `plugins` — either directly in `conf/params.config`, or (recommended) in your own `custom_config.nf`, as shown under [Customizing a run](#customizing-a-run).

## What the pipeline does

For each input VCF:

1. **`NO_G_VCF`** — strips genotypes (annotation doesn't need them, and it keeps VEP's input small).
2. **Normalisation** — `bcftools norm -m-` splits multiallelic records; if `left_align = true`, this also left-aligns indels against `ref_fasta`.
3. **Splitting** *(optional, `split_input = true`)* — splits the normalised VCF into `number_of_chunks` pieces so VEP can annotate them in parallel.
4. **`RUN_VEP`** — runs VEP with `--everything`, `--fasta`, `--dir_cache`, `--dir_plugins`, and whatever's in `plugins`.
5. **CSQ extraction** — pulls the `CSQ` INFO field back out of VEP's output into TSV form.
6. Depending on which run modes are enabled (below), the CSQ TSVs are recombined, gzipped/indexed, and either re-injected into the original VCF, published as a combined TSV, or both.

## Notes

- `left_align = true` is required whenever `wxs_tsv = true`; the pipeline will fail fast with a clear error if you enable one without the other.
- `ref_fasta` must point to an existing file whenever `left_align = true`, or the pipeline fails fast at startup rather than partway through.
- `number_of_chunks` must be a positive integer greater than 1 whenever `split_input = true`.
