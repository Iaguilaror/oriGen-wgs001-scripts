/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process vcf2tsv {

    publishDir "${params.results_dir}/02-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

	tag "$vcf"

    input:
	tuple path(vcf), path(tbi)

    output:
        path "*", emit: vcf2tsv_results

    script:
    """
    bcftools view -h $vcf | grep "ID=ANN," | cut -d":" -f2 | tr -d '"> ' | tr '|' ' ' > the_ANN_header.tmp
    bcftools view -h $vcf | tail -n1 | cut -f10- | tr "\t" " " > the_samples.tmp
    echo "CHROM POS REF ALT TYPE AF AC AN" > base.tmp

    paste -d ' ' base.tmp the_ANN_header.tmp the_samples.tmp | tr " " "\t" > header.tmp

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/TYPE\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%INFO/ANN[\t%GT]\n' \
        $vcf \
    | tr '|' '\t' > body_base.tmp

    cat header.tmp body_base.tmp > ${vcf.simpleName}.tsv

    rm *.tmp

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