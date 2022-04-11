process COMBINED_OUTPUT {
      label 'process_low'

      input:
      path(txt)
      path(log)
      path(dlog)

      output:
      path("collated_species_id_results.txt"), emit: txt
      path("collated_species_id.log"), emit: log
      path("collated_download_genomes.log"), emit: dlog

      script:
      """
      """
}