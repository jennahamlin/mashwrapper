#!/usr/bin/env nextflow
/*
========================================================================================
    mashwrapper
========================================================================================
    Github : https://github.com/jennahamlin/mashwrapper
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
        if (params.email_addy) {
            String emailMessage = "Results will be emailed to: ${params.email_addy}"
    
            if (params.email_cc) {
                emailMessage += " (CC: ${params.email_cc})"
            }
    
            println(emailMessage)

            if (params.email_cc) {
                sendMail(to: params.email_addy, 
                    cc: params.email_cc,
                    subject: " ${params.email_subject}", 
                    body: msg, 
                    attach: "${params.outdir}/combinedOutput/collated_species_id_results.txt" )
            } else {
                sendMail(to: params.email_addy,
                    subject: "${params.email_subject}",
                    body: msg,
                    attach: "${params.outdir}/combinedOutput/collated_species_id_results.txt" )
            }
        } else if (params.email_cc) {
            println("Error: The email_cc flag is specified, but email_addy is required for sending results.")
            println("Please provide a valid recipient email address using the email_addy flag.")
            println("However, your results can be viewed and the folder is called: ${c_green}${params.outdir}")
       } else {
            println(""" 
            Results will not be emailed. 
            Please check your specified out directory for the results. 
            Your results folder is called: ${c_green}${params.outdir}
            """)
       }
    }  

/*
========================================================================================
    THE END
========================================================================================
*/
