process MAKE_DATABASE {
	label 'process_medium'

	conda (params.enable_conda ? "conda-forge::python=3.7.12 conda-forge::pandas=1.3.5 conda-forge::tabulate=0.8.9 bioconda::mash=2.0" : null)
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
	' https://depot.galaxyproject.org/singularity/mulled-v2-9422771e6df1a77bc63f53d9f4428f16f50bb217:78bc1e477ae739d7d2d9bdd66e4fd3074dde5974-0' :
	'quay.io/biocontainers/mulled-v2-9422771e6df1a77bc63f53d9f4428f16f50bb217:78bc1e477ae739d7d2d9bdd66e4fd3074dde5974-0' }"

	input:
	path(msh)

	output:
	path("myMashDatabase***.msh"), emit: dmsh
	path "versions.yml", emit: versions

	script:
	"""
	currentDate=`date +"%Y-%m-%d_%T"`
	
	if ls *noMash.msh &> /dev/null; then
		rm *noMash.msh;
		echo 'Removing noMash.msh files, these were temporary files';
		if ls *.msh &> /dev/null; then
			mash sketch *.msh -o myMashDatabase.\$currentDate.msh -S 42
		else
			echo "After removal there are no *.msh files available to generate Mash database. Exiting"
			exit 1
		fi
	else
		echo 'Only .msh files in directory';
		mash sketch *.msh -o myMashDatabase.\$currentDate.msh -S 42; 
	fi

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
		  mash: \$(mash --version | sed 's/Mash //g')
	END_VERSIONS
	"""
}
