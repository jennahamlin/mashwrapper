process COMBINED_OUTPUT {
      tag "combinedoutput"
      label 'process_low'

      input:
      path(txt)
      path(log)
      file(dlog)

      output:
      path("collated_species_id_results.txt"), emit: txt
      path("collated_species_id.log"), emit: log
      path("collated_download_genomes.dlog"),  emit: dlog optional true


      script:
      """
      """
}
