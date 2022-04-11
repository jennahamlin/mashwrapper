/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMashwrapper.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.organism, params.database ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input reads samplesheet not specified!' }
if (params.organism) {ch_organism = Channel.fromPath(params.organism)
                                           .splitText()
                                           .map { it.replaceFirst(/\n/,'') }} else { exit 1, 'No input file of organisms to download provided!'}
if (params.database) { ch_database = file(params.database) } else { exit 1, 'No mash sketch is included!'}
/*

========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Local: Modules
include { DOWNLOAD_GENOMES } from '../modules/local/download_genomes'
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
//include { MASH_DIST                   } from '../modules/nf-core/modules/mash/dist/main'
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
    ch_results = Channel.empty()
    ch_log = Channel.empty()
    ch_download = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run Download_Genomes
    //
    DOWNLOAD_GENOMES(
        ch_organism
    )
    ch_download = ch_download.mix(DOWNLOAD_GENOMES.out.dlog)
    //
    // MODULE: Run Species_Id
    //
    SPECIES_ID (
        ch_database, INPUT_CHECK.out.reads
    )
    ch_results = ch_results.mix(SPECIES_ID.out.txt)
    ch_log = ch_log.mix(SPECIES_ID.out.log)

    //
    // MODULE: Collate results and log into one file to send to output
    //
    COMBINED_OUTPUT (
        ch_results.unique().collectFile(name: 'collated_species_id_results.txt'), ch_log.unique().collectFile(name: 'collated_species_id.log'), ch_download.unique().collectFile(name: 'collated_download_genomes.log')

    )
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
