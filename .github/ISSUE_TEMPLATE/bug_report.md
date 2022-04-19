
```yaml
name: Bug report
description: Report something that is broken or incorrect
labels: bug
body:

  - type: markdown
    attributes:
      value: |
        Before you post this issue, please check the documentation:

        - [nf-core website: troubleshooting](https://nf-co.re/usage/..
        - [nf-core/jennahamlin-nf-core-mashwrapper pipeline document..

  - type: textarea
    id: description
    attributes:
      label: Description of the bug
      description: A clear and concise description of what the bug is.
    validations:
      required: true

  - type: textarea
    id: command_used
    attributes:
      label: Command used and terminal output
      description: Steps to reproduce the behaviour. Please paste th..
      render: console
      placeholder: |
        $ nextflow run ...

        Some output where something broke

  - type: textarea
    id: files
    attributes:
      label: Relevant files
      description: |
        Please drag and drop the relevant files here. Create a `.zip..
        Your verbose log file `.nextflow.log` is often useful _(this..

  - type: textarea
    id: system
    attributes:
      label: System information
      description: |
        * Nextflow version _(eg. 21.10.3)_
        * Hardware _(eg. HPC, Desktop, Cloud)_
        * Executor _(eg. slurm, local, awsbatch)_
        * Container engine: _(e.g. Docker, Singularity, Conda, Podma..
        * OS _(eg. CentOS Linux, macOS, Linux Mint)_
        * Version of nf-core/jennahamlin-nf-core-mashwrapper _(eg. 1..
```
---- 

Produced by [yaml2md](https://www.google.com) and [pandoc](pandoc.org)
