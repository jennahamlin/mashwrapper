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

/*
========================================================================================
    THE END
========================================================================================
*/
