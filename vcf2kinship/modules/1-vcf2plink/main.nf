/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process vcf2plink {

    publishDir "${params.results_dir}/01-vcf2plink/", mode:"symlink"

    input:
      tuple path(vcf), path(tbi)

    output:
        path "*.converted2plink*", emit: vcf2plink_results

    script:
    """
    plink2 --vcf $vcf \
      --maj-ref force \
      --vcf-require-gt \
      --make-bed \
      --out ${vcf.simpleName}.converted2plink
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