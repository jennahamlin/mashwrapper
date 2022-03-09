process SPECIES_ID {
publishDir "${params.outdir}/speciesId",
      mode: "copy"

input:
path database
tuple val(meta), path(reads)

script:
"""
gunzip -cf "${reads[0]}" > read1.fastq
gunzip -cf "${reads[1]}" > read2.fastq

run_species_id.py -d ${database} -r1 read1.fastq -r2 read2.fastq --max_dist ${params.max_dist} --min_kmer ${params.min_kmer} --num_threads ${params.num_threads}
"""
}
