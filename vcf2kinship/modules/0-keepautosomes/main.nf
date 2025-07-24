/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process keepautosomes {

    publishDir "${params.results_dir}/0-keepautosomes/", mode:"symlink"

    input:
      tuple path(vcf), path(tbi)

    output:
        path "*.vcf.gz*", emit: keepautosomes_results

    script:
    """
    bcftools view -r \
      1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
      $vcf \
      | bcftools norm -m - \
        --multi-overlaps 0 \
        -Oz -o ${vcf.simpleName}_spltimultialt_autosomes.vcf.gz \
    && tabix -p vcf ${vcf.simpleName}_spltimultialt_autosomes.vcf.gz
    """

}

/* name a flow for easy import */
workflow KEEPAUTOSOMES {

 take:
    vcf_channel

 main:

    vcf_channel | keepautosomes 

  emit:
    keepautosomes.out[0]

}
