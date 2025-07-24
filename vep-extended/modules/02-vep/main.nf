/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process vep {

    publishDir "${params.results_dir}/02-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    tag "$vcf"

    input:
      path vcf

    output:
        path "*.vcf.gz*", emit: vep_results
		path "*.html"

    script:
    """
	# --af			ANN field: AF
 	#     Add the global allele frequency (AF) from 1000 Genomes Phase 3 data for any known co-located variant to the output.
	#     For this and all --af_* flags, the frequency reported is for the input allele only, not necessarily the non-reference
	#     or derived allele.
	# --max_af		ANN field: MAX_AF, MAX_AF_POPS
 	#	  Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD.
	# --af_1kg		ANN field: AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF
	# 	  Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output.
	#     Must be used with --cache.
	# --af_gnomade	ANN field: gnomADe_AF, gnomADe_AFR_AF, gnomADe_AMR_AF, gnomADe_ASJ_AF, gnomADe_EAS_AF, gnomADe_FIN_AF, gnomADe_NFE_AF, gnomADe_OTH_AF, gnomADe_SAS_AF
	#     Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD
	#     exomes are included;
	vep \
		--input_file $vcf \
		--debug \
		--format "vcf" \
		--output_file ${vcf.simpleName}_vep.vcf.gz \
		--compress_output bgzip \
		--vcf \
		--vcf_info_field ANN \
		--force_overwrite \
		--stats_file ${vcf.simpleName}.stats.html \
		--warning_file ${vcf.simpleName}.err.txt \
		--fork 1 \
		--species "homo_sapiens" \
		--assembly GRCh38 \
		--cache \
		--offline \
		--cache_version 113 \
		--buffer_size 10000 \
		--pick \
		--canonical \
		--biotype \
		--dont_skip \
		--clin_sig_allele 0 \
		--plugin LoF,loftee_path:/home/iaguilar/bin/loftee,human_ancestor_fa:/home/iaguilar/bin/ancestral/human_ancestor_GRCh38_autosomes.fa.gz,conservation_file:/home/iaguilar/bin/phylocsf/loftee.sql \
		--dir_plugins ~/bin/loftee

	tabix -p vcf ${vcf.simpleName}_vep.vcf.gz
    """

}

/* name a flow for easy import */
workflow VEP {

 take:
    vcf_channel

 main:

    vcf_channel | vep

  emit:
    vep.out[0]

}
