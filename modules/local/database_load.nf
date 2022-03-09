process DATABASE_LOAD {
publishDir "${params.outdir}/speciesId",
      mode: "copy"

input:
path database

output:
file "${database}.msh"


"""
echo $database
"""
}
