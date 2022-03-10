process COLLATED_RESULTS {
      label 'process_low'

      input:
      path(txt)
      //path(log)

      output:
      path("collated_results.txt")  , emit: txt
      //path("log.txt") , emit: log

      script:
      """
      """
}

//makes folder but the file is empty
//process COLLATED_RESULTS {
  //    label 'process_low'
//
  //    input:
    //  path(txt)
      //path(log)
//
  //    output:
    //  path("collated_results.txt")  , emit: txt
      //path("log.txt") , emit: log
//
  //    script:
    //  """
      //< $txt cat > collated_results.txt
    //  """
//
//}
