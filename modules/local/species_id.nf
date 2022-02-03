params:options = [:]

process SPECIES_ID {

  publishDir "${params.outdir}/speciesId",
        mode: "copy"

  input:
  path database

  output:
  file "out.txt"


  script:
  """
  run_species_id.py -d ${database} > out.txt
  """

}
