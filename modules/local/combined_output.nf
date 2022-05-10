process COMBINED_OUTPUT {
      tag "combinedoutput"
      label 'process_low'

      input:
      path(txt)
      path(log)
      path(dlog)
      file(name)

      output:
      path("collated_species_id_results.txt"), emit: txt
      path("collated_species_id.log"), emit: log
      path("collated_download_genomes.log"),  emit: dlog, optional: true
      path("database.name", includeInputs: true), emit: name, optional: true


      script:
      """
      """
}

