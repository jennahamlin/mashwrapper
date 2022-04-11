process COMBINED_OUTPUT {
      label 'process_low'

      input:
      path(txt)
      path(log)

      output:
      path("collated_results_species_id_.txt")  , emit: txt
      path("collated_species_id.log") , emit: log

      script:
      """
      """
}
