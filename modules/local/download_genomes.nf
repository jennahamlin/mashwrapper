process DOWNLOAD_GENOMES {
      label 'process_low'

      conda (params.enable_conda ? "bioconda::p7zip=15.09 conda-forge::ncbi-datasets-cli=12.20.1" : null)
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/ncbi-datasets-cli:12.20.1' :
      'quay.io/biocontainers/ncbi-datasets-cli:12.11.0' }"

      input:
      val(organism)
      val(conda)
      val(assembly)

      output:
      path("*.fna"), emit: fna
      path(".command.log"), emit: dlog
      path("versions.yml"), emit: versions

      script:
      """
      ## Original downloadGenome.sh script allows user to specify using conda
      ## environment, incorporate a boolean for conda to be used with -c flag

      ## export nf variable to tell downloadGenome script if nextflow is in use
      nf="This is script is running via NextFlow"
      export nf
      echo $assembly

      if [[ "$assembly" != false ]]; then
          ${projectDir}/bin/downloadGenome.sh -c "${conda}" -s "${organism}" -a "${assembly}"
      else
          ${projectDir}/bin/downloadGenome.sh -c "${conda}" -s "${organism}"
      fi

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          datasets: \$(datasets version | sed 's/Datasets //g')
          dataformat: \$(dataformat version | sed 's/Dataformat //g')
      END_VERSIONS
      """
}
