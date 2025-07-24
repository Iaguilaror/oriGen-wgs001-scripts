/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process dbsnp_summ {

    publishDir "${params.results_dir}/05-a-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

	tag "$vcf"

    input:
	path(vcf)

    output:
        path "*", emit: dbsnp_summ_results

    script:
    """
    # get the chr in turn
    thechr=\$(bcftools view -H ${vcf} | head -n1 | cut -f1)

    bcftools view -H \
        -v snps $vcf | cut -f3 > snps.tmp

    bcftools view -H \
        -v indels $vcf | cut -f3 > indels.tmp

    total_snp=\$(wc -l snps.tmp | cut -d " " -f1)
    novel_snp=\$(grep -c "^\\." snps.tmp)
    known_snp=\$(grep -v -c "^\\." snps.tmp)

    total_indel=\$(wc -l indels.tmp | cut -d " " -f1)
    novel_indel=\$(grep -c "^\\." indels.tmp)
    known_indel=\$(grep -v -c "^\\." indels.tmp)

    echo "chromosome dbsnp n_variants type
    \$thechr total \$total_snp snp
    \$thechr novel \$novel_snp snp
    \$thechr known \$known_snp snp
    \$thechr total \$total_indel indel
    \$thechr novel \$novel_indel indel
    \$thechr known \$known_indel indel" > ${vcf.simpleName}.dbsnp_summary.txt

    rm *.tmp
    """

}

/* name a flow for easy import */
workflow DBSNP_SUMM {

 take:
    vcf_channel

 main:

    vcf_channel | dbsnp_summ

  emit:
    dbsnp_summ.out[0]

}
