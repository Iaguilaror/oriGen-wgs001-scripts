/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process dbsnp {

    publishDir "${params.results_dir}/01-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    tag "$vcf"

    input:
      tuple path(vcf), path(tbi)

    output:
        path "*", emit: dbsnp_results

    script:
    """
    # get the chr in turn
    thechr=\$(bcftools view -H ${vcf} | head -n1 | cut -f1)

    bcftools view -O z ${params.dbsnp_ref} \$thechr > small_reference.vcf.gz && tabix -p vcf small_reference.vcf.gz

    bcftools annotate \
    -a small_reference.vcf.gz \
    -c ID \
    -o ${vcf.simpleName}_dbSNP.vcf.gz \
    -O z \
    $vcf

    rm small_reference.vcf.gz*
    """

}

/* name a flow for easy import */
workflow DBSNP {

 take:
    vcf_channel

 main:

    vcf_channel | dbsnp

  emit:
    dbsnp.out[0]

}
