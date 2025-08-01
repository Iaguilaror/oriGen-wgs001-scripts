
/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process getcov {

    publishDir "${params.results_dir}/01-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
      tuple path( cram ), path( crai ), path( bedfile )

    output:
        path "*.bed.gz", emit: getcov_results

    script:
    """
    samtools depth \
      --threads 2 \
      -b $bedfile \
      $cram \
      | sort -n -k1 -k2 - \
      | awk 'BEGIN { FS = OFS = "\t" } { print \$1, \$2, \$2+1, \$3 }' \
      > ${cram.simpleName}.readsperposition.bed \
      && bgzip --compress-level 4 \
	--threads 2 \
	${cram.simpleName}.readsperposition.bed
    """

}

/* name a flow for easy import */
workflow GETCOV {

 take:
    data_channel

 main:

    data_channel | getcov 

  emit:
    getcov.out[0]

}
