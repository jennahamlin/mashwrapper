process SPECIES_ID {

  publishDir "${params.outdir}/speciesId",
        mode: "copy"

  input:
  path database

  output:
  file "out.txt"


  script:
  """
  run_species_id.py -d ${database} --max_dist ${params.max_dist} --min_kmer ${params.min_kmer}  > out.txt
  """

}
