process MAKE_MASH {
      label 'process_low'

      //container = "file://mashpython_v1.sif"
      conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      ' https://depot.galaxyproject.org/singularity/mulled-v2-9422771e6df1a77bc63f53d9f4428f16f50bb217:78bc1e477ae739d7d2d9bdd66e4fd3074dde5974-0' :
      'quay.io/biocontainers/mulled-v2-9422771e6df1a77bc63f53d9f4428f16f50bb217:78bc1e477ae739d7d2d9bdd66e4fd3074dde5974-0' }"

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
