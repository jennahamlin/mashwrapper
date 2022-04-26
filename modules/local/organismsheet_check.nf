process ORGANISMSHEET_CHECK {
    tag "$organismsheet"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.7.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.7.6' :
        'quay.io/biocontainers/python:3.7.6' }"

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
