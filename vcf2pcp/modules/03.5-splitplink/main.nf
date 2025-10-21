/* Inititate DSL2 */
nextflow.enable.dsl=2

/* Define the main processes */
process splitplink {

    publishDir "${params.results_dir}/03.5-${task.process.replaceAll(/.*:/, '')}/", mode:"symlink"

    input:
        tuple path(one), path(two), path(three), path(four), path(reftsv), path(popsprojection)

    output:
        path "*", emit: splitplink_results

    script:
    """
    # Store the formatted string in a variable
    terms=\$(cat $popsprojection | tr "\n" "|" | sed 's/|\$//')

    # Use the variable with grep
    grep -w -E "\$terms" $reftsv | cut -f1 | awk ' BEGIN{FS=OFS="\t"} {print "ori-SAMPLE", \$1}' > REF_IDS.tmp

    grep -v -w -E "\$terms" $reftsv | grep -v "sample" | grep -v "#" | cut -f1 | awk ' BEGIN{FS=OFS="\t"} {print "ori-SAMPLE", \$1}' > NONREF_IDS.tmp

    # Input file
    input_file="NONREF_IDS.tmp"

    # Shuffle the input file and split it
    shuf "\$input_file" \
    | tee >(head -n ${params.random_samples} > "random_lines.tmp") \
    | tail -n +\$(( ${params.random_samples} + 1 )) > "rest_of_lines.tmp"

    # put the randoms in the REFs
    cat random_lines.tmp >> REF_IDS.tmp

    # Rename the remaining nonrefs
    mv rest_of_lines.tmp NONREF_IDS.tmp

    plink2 --bfile allwgs.converted2plink \
       --keep REF_IDS.tmp \
       --make-bed \
       --out references_for_admixture

    plink2 --bfile allwgs.converted2plink \
       --keep NONREF_IDS.tmp \
       --make-bed \
       --out problem_for_admixture
    """

}

/* name a flow for easy import */
workflow SPLITPLINK {

 take:
    data_channel

 main:
    data_channel | splitplink

  emit:
    splitplink.out[0]

}
