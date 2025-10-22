/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process vcf2plink {

    publishDir "${params.results_dir}/03-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    tag "$vcf"

    input:
      tuple path(vcf), path(tbi)

    output:
        path "*.converted2plink*", emit: vcf2plink_results
        path "number_of_samples_and_snps"

    script:
    """
    plink2 --vcf ${vcf} \
      --maj-ref force \
      --const-fid "ori-SAMPLE" \
      --vcf-require-gt \
      --make-bed \
      --vcf-filter \
      --vcf-half-call m \
      --out allwgs.converted2plink

    echo "# Dataset: ${vcf.simpleName}.converted2plink" > number_of_samples_and_snps.txt
    echo "n_samples n_snps" >> number_of_samples_and_snps
    echo "\$(wc -l allwgs.converted2plink.fam) \$(wc -l allwgs.converted2plink.bim)" >> number_of_samples_and_snps

    """

}

/* name a flow for easy import */
workflow VCF2PLINK {

 take:
    vcf_channel

 main:

    vcf_channel | vcf2plink 

  emit:
    vcf2plink.out[0]

}