/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process mcps {

    publishDir "${params.results_dir}/04-${task.process.replaceAll(/.*:/, '')}/", mode:"copy"

    tag "$vcf"

    input:
      tuple path(vcf), path(tbi)

    output:
        path "*", emit: mcps_results

    script:
    """
	# get the chr in turn
    thechr=\$(bcftools view -H ${vcf} | head -n1 | cut -f1)

    bcftools view -O z ${params.mcps_ref} \$thechr > small_reference.vcf.gz && tabix -p vcf small_reference.vcf.gz
	
    bcftools annotate \
    -a small_reference.vcf.gz \
    --columns INFO \
    -o ${vcf.simpleName}_MCPS.vcf.gz \
    -O z \
    $vcf

    rm small_reference.vcf.gz*
    """

}

/* name a flow for easy import */
workflow MCPS {

 take:
    vcf_channel

 main:

    vcf_channel | mcps

  emit:
    mcps.out[0]

}
