process DOWNLOAD_GENOMES {
      label 'process_low'
      
      conda (params.enable_conda ? "conda-forge::ncbi-datasets-cli=12.20.1" : null)
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/ncbi-datasets-cli:12.20.1' :
      'quay.io/biocontainers/ncbi-datasets-cli:12.11.0' }"

      input:
      val(organism)

      output:
      path("**/allDownload/*.fna"), emit: fna
      path(".command.log"), emit: dlog
      path("versions.yml"), emit: versions

      script:
      """
      ## Original downloadGenome.sh script allows user to specify using conda
      ## environment, so always set -c flag to False in Nextflow
      
      ${projectDir}/bin/downloadGenome.sh -c F -s "${organism}"

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          datasets: \$(datasets version | sed 's/Datasets //g')
          dataformat: \$(dataformat version | sed 's/Dataformat //g')
      END_VERSIONS
      """
}
