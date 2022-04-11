process DOWNLOAD_GENOMES {
      label 'process_low'

      input:
      val organism

      output:
      // is many *fna files that will be used to generate .msh and a database of .msh
      path(".command.log"), emit: dlog
      path("*.fna"), emit: fasta

      script:
      """
      ## Original downloadGenome.sh script allows user to specify using conda
      ## environment. Because of that, we set the -c always to False in Nextflow
      ${projectDir}/bin/downloadGenome.sh -c F -s "$organism"
      """
}
