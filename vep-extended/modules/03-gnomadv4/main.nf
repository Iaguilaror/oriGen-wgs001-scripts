/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process gnomad {

    publishDir "${params.results_dir}/03-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    tag "$vcf"

    input:
      tuple path(vcf), path(tbi)

    output:
        path "*", emit: gnomad_results

    script:
    """
	# get the chr in turn
    thechr=\$(bcftools view -H ${vcf} | head -n1 | cut -f1)

    bcftools view -O z ${params.gnomad_ref} \$thechr > small_reference.vcf.gz && tabix -p vcf small_reference.vcf.gz
	
    bcftools annotate \
    -a small_reference.vcf.gz \
    --columns INFO \
    -o ${vcf.simpleName}_gnomADv4.vcf.gz \
    -O z \
    $vcf

    tabix -p vcf ${vcf.simpleName}_gnomADv4.vcf.gz

    rm small_reference.vcf.gz*
    """

}

/* name a flow for easy import */
workflow GNOMAD {

 take:
    vcf_channel

 main:

    vcf_channel | gnomad

  emit:
    gnomad.out[0]

}
