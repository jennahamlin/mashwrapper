#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/mashwrapper
========================================================================================
    Github : https://github.com/nf-core/mashwrapper
    Website: https://nf-co.re/mashwrapper
    Slack  : https://nfcore.slack.com/channels/mashwrapper
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    DATABASE PARAMETER VALUES
========================================================================================
*/

//params.database = WorkflowMain.getGenomeAttribute(params, 'msh')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { MASHWRAPPER } from './workflows/mashwrapper'

//
// WORKFLOW: Run main nf-core/mashwrapper analysis pipeline
//
workflow NFCORE_MASHWRAPPER {
    MASHWRAPPER ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_MASHWRAPPER ()
}

 workflow.onComplete {
    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
        .stripIndent()

        sendMail(to: 'ptx4@cdc.gov', subject: "Results from mashWrapper for identifying species - ${params.email_subject}", body: msg, attach: "${params.outdir}/combinedOutput/collated_species_id_results.txt" )
    }

/*
========================================================================================
    THE END
========================================================================================
*/
