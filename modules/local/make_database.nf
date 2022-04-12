process MAKE_DATABASE {
      label 'process_medium'

      input:
      //file(msh).collect()
      path(msh)

      output:
      path("myMashDatabase.msh"), emit: dmsh


      script:
      """
      mash sketch $msh -o myMashDatabase.msh
      """
}
