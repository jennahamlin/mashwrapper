process MAKE_DATABASE {
      label 'process_medium'

      container = "file://mashpython_v1.sif"

      input:
      path(msh)

      output:
      path("myMashDatabase.msh"), emit: dmsh

      script:
      """
      mash sketch $msh -o myMashDatabase.msh
      """
}
