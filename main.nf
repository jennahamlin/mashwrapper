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

c_green = "\033[0;32m";
 workflow.onComplete {
    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at : ${workflow.complete}
        Duration : ${workflow.duration}
        Success : ${workflow.success}
        Command line : ${workflow.commandLine}
        workDir : ${workflow.workDir}
        Exit status : ${workflow.exitStatus}
        Error Message : ${workflow.errorMessage ?: 'None'}
        Error report : ${workflow.errorReport ?: '-'}
        Nextflow version : ${workflow.nextflow.version}
        Nextflow build : ${workflow.nextflow.build}
        Nextflow Compile Timestamp : ${workflow.nextflow.timestamp}
        """
        .stripIndent()
        if (params.email_addy && params.email_cc) {
        sendMail(to: params.email_addy, 
                 cc: params.email_cc,
                 subject: "Results from mashWrapper for identifying species - ${params.email_subject}", 
                 body: msg, 
                 attach: "${params.outdir}/combinedOutput/collated_species_id_results.txt" )
        } else { 
               if (params.email_addy && !params.email_cc) {
               sendMail(to: params.email_addy,
                        subject: "Results from mashWrapper for identifying species - ${params.email_subject}",
                        body: msg,
                        attach: "${params.outdir}/combinedOutput/collated_species_id_results.txt" )
       } else {
                println("""
                Results will not be emailed. 
                Please check your specified out directory for the results. 
                Your results folder is called: ${c_green}${params.outdir}
                """)
       }
    }
  }

/*
========================================================================================
    THE END
========================================================================================
*/
