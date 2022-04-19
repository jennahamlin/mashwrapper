process COMBINED_OUTPUT {
      label 'process_low'

      input:
      path(txt)
      path(log)
      path(dlog) 

      output:
      path("collated_species_id_results.txt"), emit: txt
      path("collated_species_id.log"), emit: log
      path("collated_download_genomes.log"),  emit: dlog optional true
      path "versions.yml", emit: versions

      script:
      """
      
      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$(python --version | sed 's/Python //g')
      END_VERSIONS
      """
}
