process MAKE_MASH {
      label 'process_high_memory'

      input:
      file(fna).collect()

      output:
      path("*.msh"), emit: msh


      script:
      """
      echo ${fna[0]}

      for file in ${fna}
      do
          mash sketch \$file -k 25 -s 100000
      done
      """
}
