process COLLATED_RESULTS {
      label 'process_low'

      input:
      path(txt)
      //path(checks)

      output:
      path("collated_results.txt")  , emit: txt
      //path("collated_checks.txt") , emit: log
      
      script:
      """
      """
}
