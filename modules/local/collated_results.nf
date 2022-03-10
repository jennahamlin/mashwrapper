process COLLATED_RESULTS {
      label 'process_low'

      input:
      path(txt)
      path(log)

      output:
      path("*.txt")  , emit: results
      path("run.log") , emit: log

      script:
      def args = task.ext.args ?: ''

}
