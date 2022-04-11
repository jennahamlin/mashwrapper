process MAKE_MASH {
      label 'process_low'

      input:
      file fasta

      output:
      // is many *fna files that will be used to generate .msh and a database of .msh
      //path(".command.log"), emit: dlog
      //path("*.fna"), emit: fasta

      script:
      """
      mash sketch "$fasta" -k 25 -s 100000 -p 8
      """
}
