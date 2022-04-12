process DOWNLOAD_GENOMES {
      label 'process_low'

      input:
      val(organism)

      output:
      // is many *fna files that will be used to generate .msh files and then mash database
      path("**/allDownload/*.fna"), emit: fna
      path(".command.log"), emit: dlog

      script:
      """
      ## Original downloadGenome.sh script allows user to specify using conda
      ## environment, so always set -c flag to False in Nextflow
      ${projectDir}/bin/downloadGenome.sh -c F -s "${organism}"
      """
}
