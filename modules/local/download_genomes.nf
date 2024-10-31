process DOWNLOAD_GENOMES {
      label 'process_low'

      conda (params.enable_conda ? "bioconda::p7zip=15.09 conda-forge::ncbi-datasets-cli=15.2.0" : null)
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'quay.io/staphb/ncbi-datasets:15.2.0' :
      'quay.io/staphb/ncbi-datasets:15.2.0' }"

      input:
      val(organism)
      val(conda)
      val(assembly)

      output:
      path("*.fna"), emit: fna
      path(".command.log"), emit: dlog
      path("versions.yml"), emit: versions
      path("excluded_genomes.txt"), emit: incon

      script:
      """
      # Export variable to indicate script is running via NextFlow
      export nf="This script is running via NextFlow"

      # Execute the downloadGenome command
      if [[ "$assembly" != false ]]; then
          ${projectDir}/bin/downloadGenome.sh -c "${conda}" -s "${organism}" -a "${assembly}"
      else
          ${projectDir}/bin/downloadGenome.sh -c "${conda}" -s "${organism}"
      fi

      # Create versions.yml file with dataset and dataformat versions
      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          datasets: \$(datasets version | awk '{print \$3}')
          dataformat: \$(dataformat version | awk '{print \$2}')
      END_VERSIONS
      """
}
