process MAKE_MASH {
      label 'process_low'

      container = "file://mashpython_v1.sif"

      input:
      file(fna).collect()

      output:
      path("*.msh"), emit: msh


      script:
      """
      for file in ${fna}
      do
          mash sketch \$file -k 25 -s 100000
      done
      """
}
