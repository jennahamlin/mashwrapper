process DOWNLOAD_GENOMES {
      label 'process_low'

      input:
      //is a list of species, need to be in quotes
      //is either true or false but I think should always be false...
      val organism


      output:
      // is many *fna files that will be used to generate .msh and a database of .msh
      //path("collated_results.txt"), emit: txt
      //path 'database'

      script:
      """
      #how to deal with if multiple -s flags are wanted.... like below
      ${projectDir}/bin/downloadGenome.sh -c F -s "$organism"
      """
}
