#!/bin/bash 

## Author J. Hamlin
## Date 01/05/2022
## This script is to download all available *legionella* and
## *legionella*-associated genomes from NCBI using the NCBI datasets command
## line tools. The NCBI dataformat command line tool then summarizes the output
## by species count. Additonally, it converts non-species descriptive file names
## (e.g., GCA_XXXX) to include genus and species in the file name. In theory,
## this tool should download data for any species located on NCBI by specifying
## the organism of interest through the -s flag. Though I have only tested the
## tool for the species listed below so careful double checking of download and
## file conversion should be done when used for other species.
## *Species Tested:*
## "legionella" "fluoribacter" 

timestamp=$(date +%d-%m-%Y_%H-%M-%S)
basefolder=$PWD

#############################
##DESCRIPTION AND HELP MENU##
#############################

Help()
{
   ## Display Help
   echo " "
   echo "Download all genomes from NCBI using their command line \
tools (datasets/dataformat). Then convert GCA#.fna files to be named \
Genus_species_GCA#.fna"
   echo
   echo "Syntax: downloadGenome.sh [-a|c|s|h]"
   echo "Options:"
   echo " -a   Optional: Assembly-level Restrict assemblies to be \
chromosome, complete_genome, contig, scaffold"
   echo " -c   Required: Activate the conda enviroment? (T/F). \
Assumes conda environment named ncbi_datasets"
   echo " -s   Required: Download genus or species. Can use multiple -s flags. Place in quotes. \
Example: -s \"legionella\" -s \"tatlockia micdadei\" "
   echo " -h   Print this help."
   echo " "
   echo "Example command: downloadGenome.sh -c F -s \"legionella pneumophila\" -a complete_genome "
   echo " "
}

## Get the options
while getopts ":c:s:ah" option; do
   case $option in
	  h) ## Display Help by calling help function
		 Help
		 exit;;
	  a) ## Assembly-level
		 assembly=$OPTARG;;
	  c) ## Activate the conda env? (T/F)
		 conda=$OPTARG;;
	  s) ## Specifiy the organism to download. Can be multiple, seperated by a
		 ## space. Using genus species, place in quotes ("tatlockia micdadei")
		 species+=("$OPTARG");;
	  :) ## This uses the colon before the c
		 echo "You supplied a flag without an argument. Exiting." >&2
		 exit 1 ;;
	  \?) ## Handle invalid options
		 echo 'Invalid Option. Please use the -h flag to display the help menu. Exiting.' >&2 
		 exit 1 ;;
   esac
done

###################
##PARSE THE FLAGS##
###################

## Make sure both -c and -s flags are included and not empty
parse() {
if [[ $conda == "" && $species == "" ]]
then
	echo 'Please supply both -c and -s flags along with their required request.
Please look at the help menu via -h, if you need additional information. Exiting.'
	exit 1
elif [[ $conda == "" ]]
then
	echo 'You did not provide the -c flag (T/F).
Please look at the help menu, if you need additional information. Exiting.'
	exit 1
elif [[ $species == "" ]]
then
	echo 'You did not provide the -s flag. Exiting.
Please look at the help menu, if you need additional information. '
	exit 1
fi
}

## Call the parse function for dealing with invalid flags
parse

## Determine what to do based on input from -c (conda) flag
## First check if conda is specifed as True
condaOrNot(){
if [[ $conda == @(True|true|T|t) && $species ]]
then
  ## Determine if conda == T request is from nextflow or not, by that I mean
  ## is the conda environment loaded  locally or requested to be used via
  ## nextflow. $nf exported from the downloadGeome.nf module, if from nextflow.
	if [[ -z "$nf" ]]
	then
		echo $nf
		echo 'The variable check to determine if next flow is running is empty.
Activating your local conda environment. Assumes conda environment is called ncbi_datasets. Continuing...'
		## I used to be able to do just the eval statement, now it does not work
		## and i don't like that this sources from somewhere that a user might not have (miniconda)...
		#eval "$(conda shell.bash hook)" ## initialize the shell to use conda
		source ~/miniconda3/etc/profile.d/conda.sh 
		conda activate ncbi_datasets
		condaAct=`echo $CONDA_DEFAULT_ENV`
 i       echo "This is your local conda enviroment that is activated:" $condaAct

		## If conda == T is not from nextflow; then confirm name of environment
		  if [[ $condaAct == 'ncbi_datasets' ]]
		  then
			echo 'These are the tools in your local conda enviroment...'
			conda list -n ncbi_datasets
		  else
			echo "This tool assumes the ncbi datasets cli conda environment is
called ncbi_datasets. Exiting."
			exit 1
		  fi
		else
		  echo $nf "and the user has requested conda. Continuing..."
	fi
elif [[ $conda == @(False|false|F|f) && $species ]]
then
  echo "Confirming both NCBI datasets/dataformat tools are available..."

  ## Check that both tools are available. If not then exit. 
  ## Checking for mash will fail because we have not gotten to that module yet 
  ## meaning that we have not downloaded the singularity image. 
  command -v dataformat >/dev/null 2>&1 || { echo >&2 "NCBI dataformat is not installed.  Exiting."; exit 1; }
  command -v datasets >/dev/null 2>&1 || { echo >&2 "NCBI datasets is not installed.  Exiting."; exit 1; }

	  ## Inform user that the tools are accessible for when conda is false
  echo "Great tools available to access NCBI and run Mash. Continuing.."
else
  echo 'Unable to activate a conda environment, find the tools in a bin folder,
or confirm that the tool is using a container. Exiting.'
  exit 1
fi
}
condaOrNot

#################
##BEGIN PROCESS##
#################
echo " "
echo "Beginning the process..."

mkdir genomesDownloaded_$timestamp
cd genomesDownloaded_$timestamp
subfolder=$basefolder/"genomesDownloaded_$timestamp"

echo "Checking your directory..."

#####################
##CHECK FILE/FOLDER##
#####################
fileCheck(){
FILE=downloadedData.tsv
if [[ -f "$FILE" ]]; then
	echo "$FILE exists, this script will overwrite to previous $FILE \
and you will lose that summary file. Exiting."
	exit 1
else
	echo "Good, a $FILE summary file does not already exist. Continuing..."
fi

FILE2=speciesCount.txt
if [[ -f "$FILE2" ]]; then
	echo "$FILE2 exists, this script will overwrite to previous $FILE2 \
and you will lose that summary file. Exiting"
	exit 1
else
	echo "Good, the $FILE2 summary file doesn't already exist. Continuing..."
fi

DIR=allDownload
if [[ -d "$DIR" ]]; then
	echo "$DIR directory exists! Please rename or remove the $DIR directory. \
Exiting."
	exit 1
else
	echo "$DIR directory does not exist, making it now and downloading will \
begin..."
	mkdir allDownload
fi
}
fileCheck

####################
##DOWNLOAD GENOMES##
####################
## Loop through species array and download genomes.
## Can change source to genbank (GCA) or refseq (GCF.)
## I use genebank option as this has a larger number of genomes.
## But with some internal checks to not included questionable genomes
## like those with an inconclusive taxonomy status. 

# echo "If assembly-level was specified, then this was the level of restriction \
# for genomes to download: $assembly"

error_handler()
{
  echo "No $assembly" files available. Creating a file place holder for this species: $valUp. Exiting.
  cd ..
  echo "There are no $assembly files available at the level specified. Exiting." > $valUp-$assembly-noFNA.fna
  echo "Exiting from this isolate..."
  echo "----------------------------------------"
}

for val in "${species[@]}";
do
  echo "This is one of the species that will be downloaded to make the mash database: ${species[@]}"
  valUp="${val:1:-1}"                                                           ## Remove quotes
  valUp=${val//[[:blank:]]/}                                                    ## Remove space

  echo "Beginning to dowload genomes from NCBI..."

	  if [[ -z "$assembly" ]] ; then
		echo "Assembly level is not specified as the parameter is empty..."

		## Changed complete_genome to complete and no longer required to specify which files to exclude
		## with datasets version change 12.2.0 -> 15.2.0
		datasets download genome taxon "$val" --dehydrated --assembly-source genbank \
--filename $valUp.zip --assembly-level complete,chromosome,scaffold,contig 

	  elif [[ -n "$assembly"  ]]; then
		echo "Assembly level is specified and will only download $assembly..."
		datasets download genome taxon "$val" --dehydrated --assembly-source genbank \
 --filename $valUp.zip --assembly-level "$assembly"   2>/dev/null || error_handler
	  fi
	   

	  ## Confirm zip file avaiable
	  count=`ls -1 *.zip 2>/dev/null | wc -l`
	  if [ i$count != 0 ]; then
		
		echo "Checking to see if there are any zip files, which there should be after downloading data..."
		echo "The number of zipped files is:"
		echo $count
		echo "Unzipping the associated files..."
		unzip $valUp.zip -d $valUp
		echo " "

		datasets rehydrate --directory $valUp

		## Checking for genomes where NCBI taxonomy is NOT OK, adding the genebank/refseq id to a list (excluded_genomes)
		## Ultimately if genome GCA ID added to list, it will be deleted 
		cat $basefolder/genomesDownloaded_$timestamp/$valUp/ncbi_dataset/data/assembly_data_report.jsonl | awk '{if (!/OK/) print $1}' | grep -o "GCA_..........." >> excluded_genomes.tmp
		
		## Get 'unculture' legionella species GCA ids and adde to a list (excluded_genomes)
		cat $basefolder/genomesDownloaded_$timestamp/$valUp/ncbi_dataset/data/assembly_data_report.jsonl |  awk '{if(/uncultured/) print $1}' | grep -o "GCA_..........." >> excluded_genomes.tmp
		
		## I want to exclude genomes which have a lower level of completeness (< 93.00), 
		## but all fields in the assembly file are not consistent between isolates. For example 
		## some do not have a contamiation estimate, so I get the data which does and then parse 
		## it here (print $1 $3) and below I get the data that does not and parse it (print $1 $4)
		## These genome GCA IDs get added to the excluded_genomes list
		## Value of 93.00 was chosen based on most isolates were higher than this value. 
		## TODO; Determine number of isolates with completeness estimate of 93.00 or higher
		cat $basefolder/genomesDownloaded_$timestamp/$valUp/ncbi_dataset/data/assembly_data_report.jsonl |  grep -o "completeness.*" | grep -o ".*organism" | awk -F , '{ print $1 $3 }' | grep -v "contamination" | awk -F\" '{ print $2 " " $5 }' | awk -F : '{ print $2 }' | awk '{ if( $1 < 93.00) print $2 }' >> excluded_genomes.tmp
 
		cat $basefolder/genomesDownloaded_$timestamp/$valUp/ncbi_dataset/data/assembly_data_report.jsonl | grep -o "completeness.*" | grep -o ".*organism" | awk -F , '{ print $1 $4 }' | grep -v "contamination" | awk -F\" '{ print $2 " " $5 }' | awk -F : '{ print $2 }' | awk '{ if( $1 < 93.00) print $1 " " $2 }' | grep GCA | awk '{ print $2 }' >> excluded_genomes.tmp

		cat excluded_genomes.tmp | uniq -u >> excluded_genomes.txt
 
		## Now remove folders for genomes listed in exclusion file
		TO_BE_DEL="excluded_genomes.txt"
		while read -r file ; do
 
		  rm -r $valUp/ncbi_dataset/data/"$file" >/dev/null 
		done < "$TO_BE_DEL"

		cp excluded_genomes.txt $basefolder
		rm excluded_genomes.tmp
		
		cd $valUp/ncbi_dataset/data
		  
		## Check for plasmids and remove
		## Get fasta header /^>/ 
		## Assign header with plasmid to p
		## Remove plasmid from genome by specifying not P (!p)
		echo "Checking for plasmids. This can take a some time..."
		for i in */*.fna
		do
		  awk '/^>/ { p = ($0 ~ /plasmid/) } !p' $i > ${i%\.*}_cleaned.fna
		  mv ${i%\.*}_cleaned.fna $i
		done

		find . -name "chrunnamed*.unlocalized.scaf.fna" -exec rm -rf {} \;				## These are plasmid files also
		find . -name "*.fna" -exec grep "plasmid" {} \; -exec rm {} \;
		find . -name "cds_from_genomic.fna" -exec rm -rf {} \;							## Remove these files. Downloaded via conda
		find . -size 0 -type f -delete													## Remove files with zero bytes

		## Make summary file of the downloaded data
		## Get opposite of grep, so do not pass excluded genome information for file manipulation/checking
		echo "Making $valUp map file for file name conversion and converting file names..."

		dataformat tsv genome --inputfile *.jsonl --fields organism-name,accession,assminfo-paired-assmaccession --force  | \
		grep -vFwf $basefolder/genomesDownloaded_$timestamp/excluded_genomes.txt  | \
		awk 'FNR==1 { header = $0; print }  $0 != header' | \
		sed 's/\//-/g' | \
		sed 's/ /_/g' | tee temp | \
		awk '{ print $1 "_" $2 }' > temp2 				## Now combine column 1 with underscore and column 2

		## Remove headers
		sed -i '1d' temp ##downloaded-$valUp.tsv
		sed -i '1d' temp2 ##map2$valUp.txt

		## Make final map file
		cut -f2 temp | paste -d " " temp2 - | \
		awk '{ print $2 ".fna" " " $1 }' > mapFinal$valUp.txt 				## Makes file of two columns old file name and new file name
        
		## Rename *.fna to folderName.fna to handle unplaced scaffolds and duplicate file names between isolates
		for d in */
		do
		  FILEPATH=$d*.fna 
		  mv $FILEPATH "$(dirname "$FILEPATH")/$(dirname "$FILEPATH").fna"
		done
        
		cp */*.fna $basefolder/genomesDownloaded_$timestamp/allDownload
		
		## Rename files using mapfile
		cd $basefolder/genomesDownloaded_$timestamp/allDownload 
		awk -F " " 'system("mv " $1 " " $2 ".fna")'  $basefolder/genomesDownloaded_$timestamp/$valUp/ncbi_dataset/data/mapFinal$valUp.txt
		
		## Move back to base directory of genomesDownloaded_timestamp
		cd $subfolder
		echo "Summarizing the entire download..."

####################
##CREATE SUMMARIES##
####################
## INPUT TSV file form NCBI command line tool dataformat 
## OUTPUT TSV file with species and genebank accession downloaded
		
		dataformat tsv genome --package $valUp.zip --fields organism-name,accession,assminfo-paired-assmaccession --force | \
		grep -vFwf $basefolder/genomesDownloaded_$timestamp/excluded_genomes.txt | \
		awk 'FNR==1 { header = $0; print }  $0 != header' | \
		grep -v 'Legionella sp\. ' | \
		grep -v 'Legionella sp\. ' >> downloadedData.tsv    ## Must include space before or get rid of Lp unknown species

		## Create a txt file of a count of all species downloaded
		cat downloadedData.tsv | sed 1d | cut -f1 | cut -f2 -d ' ' | sort | uniq -c > speciesCount.txt

		## Get total number of isolates; store and compare to those in allDownload directory
		countFile=$(awk '{ SUM+=$1}END{print SUM }' speciesCount.txt)
		awk '{ SUM+=$1}END{print SUM " Total Isolates"}' speciesCount.txt >> speciesCount.txt

## TODO clean up
		cd allDownload || { echo "Failed to change directory to allDownload"; exit 1; }

		## Exclude legionella that is not identified to species and endosymbionts
		if [[ "${species^}" == "Legionella" ]]; then
		  echo "Removing Legionella only identified to genus (e.g., Legionella sp.)..."
		  rm Legionella_sp._*.fna
		else
		  echo "No extra files to remove for non-Legionella species."
		fi

		## Count number of files in folder with those in speicesCount file for comparison
		countFolder=$(ls | wc -l)

		if [[ $countFile  -eq  $countFolder ]]; then
			echo "The number of isolates in the speciesCount file matches the number of files in allDownload directory.";
		else
			echo "Mismatch between speciesCount file and downloaded files. Investigate or retry the script.";
			exit 1
		fi

		## Move files up to basefolder to all easier copying via nextflow process
		count=$(ls -1 *.fna 2>/dev/null | wc -l)
		if [[ $count -gt 0 ]]; then
			mv *.fna "$basefolder" || echo "Failed to move .fna files."    
			#rm -rf $basefolder/genomesDownloaded_$timestamp/$valUp
		else
		   echo "No .fna files generated"
		fi
		echo "Exiting the program."
		echo "-----------------------------"
  fi
done
