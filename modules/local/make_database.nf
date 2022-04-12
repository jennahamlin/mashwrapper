process MAKE_DATABASE {
      label 'process_medium'

      input:
      //file(msh).collect()
      path(msh)

      output:
      path("myMashDatabase.msh"), emit: dmsh


      script:
      """
      echo ${msh}
      mash sketch $msh -o myMashDatabase.msh
      """
}
