process MAKE_MASH {
      label 'process_low'

      conda (params.enable_conda ? "conda-forge::python=3.7.12 conda-forge::pandas=1.3.5 conda-forge::tabulate=0.8.9 bioconda::mash=2.0" : null)
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      ' https://depot.galaxyproject.org/singularity/mulled-v2-9422771e6df1a77bc63f53d9f4428f16f50bb217:78bc1e477ae739d7d2d9bdd66e4fd3074dde5974-0' :
      'quay.io/biocontainers/mulled-v2-9422771e6df1a77bc63f53d9f4428f16f50bb217:78bc1e477ae739d7d2d9bdd66e4fd3074dde5974-0' }"

      input:
      file(fna).collect()
      val(kmer)

      output:
      path("*.msh"), emit: msh
      path "versions.yml", emit: versions


      script:
      """
      for file in ${fna}
      do
      if  [[ "${fna}" != "*-noFNA.fna" ]] ; then
        mash sketch \$file -k ${kmer} -s 100000 -S 42
      else
        echo "no .fna files found for the isolate, so no mash file will be generated" > noMash_"${fna}".msh
      fi
      done

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          mash: \$(mash --version | sed 's/Mash //g')
      END_VERSIONS
      """
}
