process SPECIES_ID {
      tag "$meta.id"
      label 'process_low'

      container = "file://mashpython_v1.sif"

      input:
      path database
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

      ${projectDir}/bin/run_species_id.py -d ${database} -r1 "\${readsIn0%.gz}"  -r2 "\${readsIn1%.gz}" --max_dist ${params.max_dist} --min_kmer ${params.min_kmer} --num_threads ${params.num_threads}
      """
}
