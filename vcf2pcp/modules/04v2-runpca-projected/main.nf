/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process smartpca {

    publishDir "${params.results_dir}/04v2-projected-${task.process.replaceAll(/.*:/, '')}/", mode:"copyNoFollow"

    input:
      path MATERIALS

    output:
        path "*", emit: smartpca_results

    script:
    """
    echo "genotypename: allwgs.converted2plink.bed" > smartpca.par
	  echo "snpname: allwgs.converted2plink.bim" >> smartpca.par
	  echo "indivname: allwgs.pedind" >> smartpca.par
	  echo "evecoutname: allwgs.evec" >> smartpca.par
	  echo "evaloutname: allwgs.eval" >> smartpca.par
	  echo "altnormstyle: NO" >> smartpca.par
	  echo "numoutlieriter: 0" >> smartpca.par
	  echo "numoutevec: 10" >> smartpca.par
    echo "poplistname: refpops-for-projection.txt"  >> smartpca.par
    #
    smartpca -p smartpca.par > smartpca.stdout
    #
    grep -A 10 "^## Tracy-Widom statistics" smartpca.stdout \
	  | tr -s " " | tr " " "\t" | sed 's#\t##' > allwgs.tracy_widom_statistics
    #
    cat allwgs.evec | tr -s " " | sed 's/^ //' | tr -s " " "\t" | zip > allwgs.evec.gz \
    && rm allwgs.evec
    """

}

/* name a flow for easy import */
workflow SMARTPCA {

 take:
    data_channel

 main:

    data_channel | smartpca 

  emit:
    smartpca.out[0]

}
