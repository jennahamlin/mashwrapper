process SPECIES_ID {
      tag "$meta.id"
      label 'process_low'

      conda (params.enable_conda ? "conda-forge::python=3.7.12 conda-forge::pandas=1.3.5 conda-forge::tabulate=0.8.8 bioconda::mash=2.0" : null)
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      ' https://depot.galaxyproject.org/singularity/mulled-v2-9422771e6df1a77bc63f53d9f4428f16f50bb217:78bc1e477ae739d7d2d9bdd66e4fd3074dde5974-0' :
      'quay.io/biocontainers/mulled-v2-9422771e6df1a77bc63f53d9f4428f16f50bb217:78bc1e477ae739d7d2d9bdd66e4fd3074dde5974-0' }"

      input:
      file(inDatabase)
      tuple val(meta), path(reads)
      val(kmer)


      output:
      path("*.txt"), emit: txt
      path("*.log"), emit: log
      path("database.info"), emit: info, optional: true
      path("versions.yml"), emit: versions

      script:
      """
      kSize=\$(mash info  ${inDatabase} | awk 'FNR == 3 {print \$3}')
      export kSize

      ## converts .fastq.gz file to .fastq
      gunzip -f "${reads[0]}"
      gunzip -f "${reads[1]}"
      
      
      ## stores in variable so we can strip off .gz in the command below
      readsIn0="${reads[0]}"
      readsIn1="${reads[1]}"

      ## with work directory in /scicomp/scratch/userID, data files had no permissions
      ## changing permissions on files
      #chmod 755 "${reads[0]}"
      #chmod 755 "${reads[1]}"
      chmod 755 "\${readsIn0%.gz}"
      chmod 755 "\${readsIn1%.gz}"

      ${projectDir}/bin/run_species_id.py -b ${inDatabase} -r1 "\${readsIn0%.gz}"  -r2 "\${readsIn1%.gz}" -d ${params.max_dist} -m ${params.kmer_min} -p ${params.num_threads}

      echo $inDatabase >  "database.info"
      mash info $inDatabase >> "database.info"

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$(python --version | sed 's/Python //g')
          mash: \$(mash --version | sed 's/Mash //g')
      END_VERSIONS
      """
}
