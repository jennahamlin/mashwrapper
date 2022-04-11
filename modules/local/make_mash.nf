process MAKE_MASH {
      label 'process_low'

      input:
      file(fasta).collect()

      output:

      script:
      """
      echo ${fasta[0]}

      for file in ${fasta}
      do
          mash sketch \$file -k 25 -s 100000 -p 8
      done
      """
}
