/* Inititate DSL2 */
nextflow.enable.dsl=2

/* load functions for testing env */
// NONE

/* define the fullpath for the final location of the outs */
params.intermediates_dir = params.results_dir = "test/results"

/* load workflows for testing env */
include { KING }    from './main.nf'

/* declare input channel for testing */
// Manually pair the bed files with their other files
plink_files = Channel
        .fromFilePairs("test/data/*.{bed,fam,bim,log}", size: 4)
        .map { baseName, files -> 
            def bed = files.find { it.name.endsWith('.bed') }
            def fam = files.find { it.name.endsWith('.fam') }
            def bim = files.find { it.name.endsWith('.bim') }
            def log = files.find { it.name.endsWith('.log') }
            [bed, fam, bim, log]
        }

/* declare scripts channel for testing */
// NONE

workflow {
  KING( plink_files )
}