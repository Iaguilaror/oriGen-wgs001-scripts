/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process vcf2tsv {

    publishDir "${params.results_dir}/05-b-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

	tag "$vcf"

    input:
	path(vcf)

    output:
        path "*", emit: vcf2tsv_results

    script:
    """

    bcftools query \
        -f '%CHROM %POS %REF %ALT{0} %INFO/AF_ori %INFO/AC_ori %INFO/AN_ori %INFO/grpmax %INFO/AF_grpmax %INFO/AF_afr %INFO/AF_amr %INFO/AF_eas %INFO/AF_nfe %INFO/AF_sas %INFO/AF_RAW %INFO/AF_IMX\n' \
        $vcf \
    | awk ' BEGIN { print "CHROM POS REF ALT oriGen_AF oriGen_AC oriGen_AN gnomad_grpmax gnomad_AF_grpmax gnomad_AF_afr gnomad_AF_amr gnomad_AF_eas gnomad_AF_nfe gnomad_AF_sas mcps_AF_RAW mcps_AF_IMX"}
    {print \$0}' | tr " " "\t" > ${vcf.simpleName}.allelefrequencies.tsv
    """

}

/* name a flow for easy import */
workflow VCF2TSV {

 take:
    vcf_channel

 main:

    vcf_channel | vcf2tsv

  emit:
    vcf2tsv.out[0]

}
