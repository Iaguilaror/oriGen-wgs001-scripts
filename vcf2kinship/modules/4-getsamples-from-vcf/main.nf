/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process removesamples {

    publishDir "${params.results_dir}/04-getsamples-from-vcf/", mode:"copyNoFollow"

    input:
      tuple path(vcf), path(tbi)
      path samples_to_remove

    output:
        path "*", emit: vcf2plink_results

    script:
    """
    bcftools view \
      -S ^$samples_to_remove \
      -o ${vcf.simpleName}.only_unrelated_samples.vcf.gz \
      -O z $vcf \
    && tabix -p vcf ${vcf.simpleName}.only_unrelated_samples.vcf.gz

    """

}

/* name a flow for easy import */
workflow REMOVESAMPLES {

 take:
    vcf_channel
    samples_to_remove

 main:

    removesamples(vcf_channel, samples_to_remove )

  emit:
    removesamples.out[0]

}