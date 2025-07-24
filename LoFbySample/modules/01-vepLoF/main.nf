/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process vep {

    publishDir "${params.results_dir}/01-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    tag "$vcf"

    input:
      tuple path(vcf), path(tbi)

    output:
        path "*.vcf.gz*", emit: vep_results
        path "*.stats.html"

    script:
    """
    vep \
      --input_file $vcf \
      --debug \
      --format "vcf" \
      --output_file ${vcf.simpleName}_vep.vcf.gz \
      --compress_output bgzip \
      --vcf \
      --vcf_info_field ANN \
      --force_overwrite \
      --stats_file ${vcf.simpleName}.stats.html \
      --warning_file ${vcf.simpleName}.err.txt \
      --fork 1 \
      --species "homo_sapiens" \
      --assembly GRCh38 \
      --cache \
      --offline \
      --cache_version 113 \
      --buffer_size 10000 \
      --pick \
      --canonical \
      --biotype \
      --dont_skip \
      --clin_sig_allele 0 \
      --plugin LoF,loftee_path:${params.loftee_path},human_ancestor_fa:${params.human_ancestor_fa},conservation_file:${params.conservation_file} \
      --dir_plugins ~/bin/loftee

      tabix -p vcf ${vcf.simpleName}_vep.vcf.gz
    """

}

/* name a flow for easy import */
workflow VEP {

 take:
    vcf_channel

 main:

    vcf_channel | vep

  emit:
    vep.out[0]

}
