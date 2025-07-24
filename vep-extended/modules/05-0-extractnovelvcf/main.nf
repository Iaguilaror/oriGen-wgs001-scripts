/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process extractnovel {

    publishDir "${params.results_dir}/05-0-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

	tag "$vcf"

    input:
	path(vcf)

    output:
        path "*", emit: extractnovel_results

    script:
    """
    bcftools view --novel $vcf \
    | grep -v "AF_RAW" \
    | bgzip > ${vcf.simpleName}.novel.vcf.gz \
    && tabix -p vcf ${vcf.simpleName}.novel.vcf.gz
    """

}

/* name a flow for easy import */
workflow EXTRACTNOVEL {

 take:
    vcf_channel

 main:

    vcf_channel | extractnovel

  emit:
    extractnovel.out[0]

}
