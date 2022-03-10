process SPECIES_ID {
      tag "$meta.id"
      label 'process_low'

      conda (params.enable_conda ? "bioconda::mash=2.3" : null)
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'quay.io/biocontainers/mash:2.3--he348c14_1' }"

      input:
      path database
      tuple val(meta), path(reads)

      output:
      path("*.txt")  , emit: txt
      path("*.log") , emit: log

      script:
      """
      ## converts .fastq.gz file to .fastq
      gunzip -f "${reads[0]}"
      gunzip -f "${reads[1]}"

      ## stores in variable so we can strip off .gz in the command below
      readsIn0="${reads[0]}"
      readsIn1="${reads[1]}"

      ${projectDir}/bin/run_species_id.py -d ${database} -r1 "\${readsIn0%.gz}"  -r2 "\${readsIn1%.gz}" --max_dist ${params.max_dist} --min_kmer ${params.min_kmer} --num_threads ${params.num_threads}
      """
}
