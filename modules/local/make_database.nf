process MAKE_DATABASE {
      label 'process_high_memory'

      input:
      file(msh).collect()

      output:

      script:
      """
      echo ${msh[0]}

      mash dist $msh
      """
}
