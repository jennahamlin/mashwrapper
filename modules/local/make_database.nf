process MAKE_DATABASE {
      label 'process_low'

      input:
      file(msh).collect()

      output:

      script:
      """
      echo ${msh[0]}

      mash dist $msh 
      """
}
