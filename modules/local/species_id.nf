params:options = [:]

process SPECIES_ID {

  input:
  path database

  output:
  file "out.txt"


  script:
  """
  run_species_id.py -d ${database}
  """

}
