// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPECIES_ID {

  publishDir "${params.outdir}/speciesId",
        mode: "copy"

  input:
  path database

  output:
  file "out.txt"


  script:
  """
  run_species_id.py -d ${database} --max_dist ${params.max_dist} > out.txt
  """

}
