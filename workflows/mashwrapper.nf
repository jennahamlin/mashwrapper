/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMashwrapper.initialise(params, log)

// Check input path parameters to see if they exist  params.get_database, params.use_database
def checkPathParamList = [ params.input]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input reads samplesheet not specified!' }

// Check optional parameters
// params.organism is optional and will assume that the genomes need to get downloaded and the database built
if (params.get_database) {ch_get_database = Channel.fromPath(params.get_database)
                                           .splitText()
                                           .map { it.replaceFirst(/\n/,'') }} else { 'No input file of organisms to download provided!'}
if (params.use_database) {ch_inDatabase = file (params.use_database)}
/*

========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Local: Modules
include { DOWNLOAD_GENOMES } from '../modules/local/download_genomes'
include { MAKE_MASH } from '../modules/local/make_mash'
include { MAKE_DATABASE } from '../modules/local/make_database'
include { SPECIES_ID } from '../modules/local/species_id'
include { COMBINED_OUTPUT } from '../modules/local/combined_output'


// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
//def multiqc_report = []

workflow MASHWRAPPER {

    ch_versions = Channel.empty()
    ch_download = Channel.empty()
    ch_fna = Channel.empty()
    ch_msh = Channel.empty()
    ch_results = Channel.empty()
    ch_log = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate, and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    if (params.get_database) {

      //
      // MODULE: Run Download_Genomes
      //
      DOWNLOAD_GENOMES( ch_get_database )

      ch_download = ch_download.mix(DOWNLOAD_GENOMES.out.dlog)
      ch_fna = ch_fna.mix(DOWNLOAD_GENOMES.out.fna)

      //
      // MODULE: Make individual mash files for all genomes downloaded
      //
      MAKE_MASH( ch_fna )

      ch_msh = ch_msh.mix(MAKE_MASH.out.msh).collect()

      //
      // MODULE: Build mash database from individual mash files
      //
      MAKE_DATABASE( ch_msh )
      ch_inDatabase = MAKE_DATABASE.out.dmsh

      //
      // MODULE: Run Species_Id
      //
      SPECIES_ID (
          ch_inDatabase, INPUT_CHECK.out.reads
      )
      ch_results = ch_results.mix(SPECIES_ID.out.txt)
      ch_log = ch_log.mix(SPECIES_ID.out.log)

      //
      // MODULE: Collate results and log into one file to send to output
      //
      COMBINED_OUTPUT (
          ch_results.unique().collectFile(name: 'collated_species_id_results.txt'), ch_log.unique().collectFile(name: 'collated_species_id.log'), ch_download.unique().collectFile(name: 'collated_download_genomes.log')
          )
    } else {
      if (params.use_database) {
      //
      // MODULE: Run Species_Id
      //
      SPECIES_ID (
          ch_inDatabase, INPUT_CHECK.out.reads
      )
      ch_results = ch_results.mix(SPECIES_ID.out.txt)
      ch_log = ch_log.mix(SPECIES_ID.out.log)

      //
      // MODULE: Collate results and log into one file to send to output
      //
      COMBINED_OUTPUT (
          ch_results.unique().collectFile(name: 'collated_species_id_results.txt'), ch_log.unique().collectFile(name: 'collated_species_id.log'), ch_download.unique().collectFile(name: 'collated_download_genomes.log')
          )
      }
    }

/*
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
*/
    //
    // MODULE: MultiQC
    //
    //  workflow_summary    = WorkflowMashwrapper.paramsSummaryMultiqc(workflow, summary_params)
    //  ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect()
    //)
    //multiqc_report = MULTIQC.out.report.toList()
    //ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

//workflow.onComplete {
//    if (params.email || params.email_on_fail) {
//        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//    }
//    NfcoreTemplate.summary(workflow, params, log)
//}

/*
========================================================================================
    THE END
========================================================================================
*/
