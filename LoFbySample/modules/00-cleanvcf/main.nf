/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process cleanvcf {

    publishDir "${params.results_dir}/00-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    tag "$vcf"

    input:
      tuple path(vcf), path(tbi)

    output:
        path "*", emit: cleanvcf_results

    script:
    """
    bcftools +fill-tags $vcf -- -t TYPE \
    | bcftools annotate \
      -x ID,QUAL,FILTER,FORMAT,^INFO/AC,INFO/AN,INFO/AF,INFO/TYPE \
    | bcftools view \
      -e 'ALT="*"' \
    | bgzip --threads 2 > ${vcf.simpleName}_clean.vcf.gz \
    && tabix -p vcf ${vcf.simpleName}_clean.vcf.gz


    """

}

/* name a flow for easy import */
workflow CLEANVCF {

 take:
    vcf_channel

 main:

    vcf_channel | cleanvcf

  emit:
    cleanvcf.out[0]

}