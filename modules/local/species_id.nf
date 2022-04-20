process SPECIES_ID {
      tag "$meta.id"
      label 'process_low'

      //container = "file://mashpython_v1.sif"
      conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      ' https://depot.galaxyproject.org/singularity/mulled-v2-9422771e6df1a77bc63f53d9f4428f16f50bb217:78bc1e477ae739d7d2d9bdd66e4fd3074dde5974-0' :
      'quay.io/biocontainers/mulled-v2-9422771e6df1a77bc63f53d9f4428f16f50bb217:78bc1e477ae739d7d2d9bdd66e4fd3074dde5974-0' }"

      input:
      file(inDatabase)
      //path database
      tuple val(meta), path(reads)

      output:
      path("*.txt"), emit: txt
      path("*.log"), emit: log

      script:
      """
      ## converts .fastq.gz file to .fastq
      gunzip -f "${reads[0]}"
      gunzip -f "${reads[1]}"

      ## stores in variable so we can strip off .gz in the command below
      readsIn0="${reads[0]}"
      readsIn1="${reads[1]}"

      ${projectDir}/bin/run_species_id.py -d ${inDatabase} -r1 "\${readsIn0%.gz}"  -r2 "\${readsIn1%.gz}" --max_dist ${params.max_dist} --min_kmer ${params.min_kmer} --num_threads ${params.num_threads}
      """
}
