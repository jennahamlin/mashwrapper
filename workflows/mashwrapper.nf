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

if (params.get_database) {

  ch_get_database = Channel.fromPath(params.get_database)

  // Channel.fromPath(params.get_database)
  //                                       .splitText()
  //                                       .map { it.replaceFirst(/\n/,'') }
} else {

  if (params.use_database) {

  ch_inDatabase = file(params.use_database)

    } else {
     exit("""
      ERROR!
      You must specify a flag to either: use a database (--use_database) or generate the database (--get_database flag).

          --use_database: User provides path to a pre-built mash database. Assumes sketch size of and XXX
          --get_database: User provided text file of organisms to download written with on organism per row.
                          Can either be written as genus species (legionall pneumophila) or genus (legionella)
      """)
  }
}

/*

========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Local: Modules
include { ORGANISMSHEET_CHECK} from '../modules/local/organismsheet_check'
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
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

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
    INPUT_CHECK ( ch_input )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    if (params.get_database) {

      //
      // MODULE
      //
      ORGANISMSHEET_CHECK (ch_get_database)

      ch_organism = ORGANISMSHEET_CHECK.out.txt
                                            .splitText()
                                            .map { it.replaceFirst(/\n/,'') }
      ch_versions = ch_versions.mix(ORGANISMSHEET_CHECK.out.versions)

      //
      // MODULE: Run Download_Genomes
      //
      DOWNLOAD_GENOMES( ch_organism )

      ch_download = ch_download.mix(DOWNLOAD_GENOMES.out.dlog)
      ch_fna = ch_fna.mix(DOWNLOAD_GENOMES.out.fna)
      ch_versions = ch_versions.mix(DOWNLOAD_GENOMES.out.versions)

      //
      // MODULE: Make individual mash files for all genomes downloaded
      //
      MAKE_MASH( ch_fna )

      ch_msh = ch_msh.mix(MAKE_MASH.out.msh).collect()
      ch_versions = ch_versions.mix(MAKE_MASH.out.versions)

      //
      // MODULE: Build mash database from individual mash files
      //
      MAKE_DATABASE( ch_msh )

      ch_inDatabase = MAKE_DATABASE.out.dmsh
      ch_versions = ch_versions.mix(MAKE_DATABASE.out.versions)

      //
      // MODULE: Run Species_Id
      //
      SPECIES_ID ( ch_inDatabase, INPUT_CHECK.out.reads )

      ch_results = ch_results.mix(SPECIES_ID.out.txt)
      ch_log = ch_log.mix(SPECIES_ID.out.log)
      ch_versions = ch_versions.mix(SPECIES_ID.out.versions)

      //
      // MODULE: Collate results and log into one file to send to output
      //
      COMBINED_OUTPUT ( ch_results.unique().collectFile(name: 'collated_species_id_results.txt'), ch_log.unique().collectFile(name: 'collated_species_id.log'), ch_download.unique().collectFile(name: 'collated_download_genomes.log') )
      ch_versions = ch_versions.mix(COMBINED_OUTPUT.out.versions)
      
      CUSTOM_DUMPSOFTWAREVERSIONS ( ch_versions.unique().collectFile(name: 'collated_versions.yml' ) )

    } else {

      if (params.use_database) {

      //
      // MODULE: Run Species_Id
      //
      SPECIES_ID ( ch_inDatabase, INPUT_CHECK.out.reads )

      ch_results = ch_results.mix(SPECIES_ID.out.txt)
      ch_log = ch_log.mix(SPECIES_ID.out.log)
      ch_versions = ch_versions.mix(SPECIES_ID.out.versions)

      //
      // MODULE: Collate results and log into one file to send to output
      // TO DO: would like to save the name of the database that was provided and send that to the combinedOutput folder
      COMBINED_OUTPUT ( ch_results.unique().collectFile(name: 'collated_species_id_results.txt'), ch_log.unique().collectFile(name: 'collated_species_id.log'), ch_inDatabase )
      
      ch_versions = ch_versions.mix(COMBINED_OUTPUT.out.versions)
      
      CUSTOM_DUMPSOFTWAREVERSIONS ( ch_versions.unique().collectFile(name: 'collated_versions.yml' )
      )
    }
  }

    //
    // MODULE: MultiQC
    //
    //  workflow_summary    = WorkflowMashwrapper.paramsSummaryMultiqc(workflow, summary_params)
    //  ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
     //ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
     //ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
     //ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
     //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
     //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

     //MULTIQC (
    //     ch_multiqc_files.collect()
  //  )
  //  multiqc_report = MULTIQC.out.report.toList()
  //  ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
