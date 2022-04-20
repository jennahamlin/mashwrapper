process ORGANISMSHEET_CHECK {
    tag "$organismsheet"
    label 'process_low'

// does require python, will need to add container for. container should be consistent where possible
//    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
//        'quay.io/biocontainers/python:3.8.3' }"

    input:
    file(organismsheet)

    output:
    path '*.txt', emit: txt
    path "versions.yml", emit: versions

    """
    check_organismsheet.py \\
        "${organismsheet}"\\
        organismsheet_valid.txt \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
