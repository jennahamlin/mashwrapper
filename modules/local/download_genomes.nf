process DOWNLOAD_GENOMES {
      label 'process_low'

      //conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)
      //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      //'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.1--pyh5e36f6f_0' :
      //'quay.io/biocontainers/ncbi-genome-download:0.3.1--pyh5e36f6f_0'}"

      input:
      val(organism)

      output:
      // is many *fna files that will be used to generate .msh files and then mash database
      path("**/allDownload/*.fna"), emit: fna
      path(".command.log"), emit: dlog
      path "versions.yml", emit: versions

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
