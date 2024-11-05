#!/bin/bash 

## Author: J. Hamlin
## Date: 01/05/2022
## Script to download *legionella* and *legionella*-associated genomes from NCBI using command line tools.

timestamp=$(date +%d-%m-%Y_%H-%M-%S)
basefolder=$PWD

# Function to display help
Help() {
	echo " "
	echo "Download all genomes from NCBI and convert file names to Genus_species_GCA#.fna"
	echo "Syntax: downloadGenome.sh [-a|c|s|h]"
	echo "Options:"
	echo " -a   Optional: Assembly-level restrictions: complete, scaffold, contig, chromosome"
	echo " -c   Required: Activate the conda environment? (T/F). Assumes conda environment named ncbi_datasets."
	echo " -s   Required: Download genus or species. Can use multiple -s flags. Place in quotes."
	echo " -h   Print this help."
	echo "Example command: downloadGenome.sh -c F -s \"legionella pneumophila\" -a complete"
	echo " "
}

# Parse command-line options
while getopts ":c:s:a:h" option; do
	case $option in
		h) Help; exit;;
		a) assembly=$OPTARG;;
		c) conda=$OPTARG;;				## outside of nextflow
		s) species+=("$OPTARG");;
		:) echo "Flag without an argument. Exiting." >&2; exit 1;;
		\?) echo 'Invalid Option. Use -h for help. Exiting.' >&2; exit 1;;
	esac
done

###################
##PARSE THE FLAGS##
###################
# Validate required flags
if [[ -z $conda || -z $species ]]; then
	echo 'Please supply both -c and -s flags. Use -h for help. Exiting.'
	exit 1
fi

## Conda should be false, because running in NF
condaOrNot() {
	# Check if conda is false and species is set
	if [[ $conda == @(False|false|F|f) && $species ]]; then
		echo $nf
		echo "Confirming both NCBI datasets/dataformat tools are available..."

		# Check that both tools are available. If not, then exit. 
		command -v dataformat >/dev/null 2>&1 || { echo >&2 "NCBI dataformat is not installed. Exiting."; exit 1; }
		command -v datasets >/dev/null 2>&1 || { echo >&2 "NCBI datasets is not installed. Exiting."; exit 1; }

		# Inform user that the tools are accessible for when conda is false
		echo "Great tools available to access NCBI and run Mash. Continuing..."
		return  # Exit the function successfully

	# Check if conda is true and species is set
	elif [[ $conda == @(True|true|T|t) && $species ]]; then
		echo $nf
		echo "Conda is set to true, which asssumes is requested by Nextflow. Continuing..."
		
	else
		echo 'Unable to activate a conda environment in Nextflow; 
		check the bin folder or confirm container usage. Exiting.'
		exit 1
	fi
}
condaOrNot

#################
##BEGIN PROCESS##
#################
echo "Beginning the process..."

mkdir genomesDownloaded_$timestamp
cd genomesDownloaded_$timestamp
subfolder=$basefolder/"genomesDownloaded_$timestamp"

echo "Checking your directory..."

#####################
##CHECK FILE/FOLDER##
#####################
fileCheck() {
	local files=("downloadedData.tsv" "speciesCount.txt")
	local dir="allDownload"

	for file in "${files[@]}"; do
		if [[ -f "$file" ]]; then
			echo "$file exists. This script will overwrite it, and you will lose the summary file. Exiting."
			exit 1
		else
			echo "Good, $file does not already exist. Continuing..."
		fi
	done

	if [[ -d "$dir" ]]; then
		echo "$dir directory exists! Please rename or remove the $dir directory. Exiting."
		exit 1
	else
		echo "$dir directory does not exist. Creating it now and downloading will begin..."
		mkdir "$dir"
	fi
}
fileCheck

####################
##DOWNLOAD GENOMES##
####################
## Loop through species array to download genomes.
## The source can be changed to GenBank (GCA) or RefSeq (GCF).
## I prefer using the GenBank option due to its larger genome selection.
## Internal checks will exclude questionable genomes, 
## such as those with inconclusive taxonomy status.

error_handler_assembly()
{
	echo "No $assembly" files available. Creating a file place holder for: $species. Exiting.
	cd ..
	echo "There are no $assembly files available at the level specified. Exiting." > $valUp-$assembly-noFNA.fna
	touch "excluded_genomes.txt"
	echo "No genome data for $species at the $assembly genome level"  >> excluded_genomes.txt
	echo "Exiting from this isolate..."
	echo "----------------------------------------"
}

for val in "${species[@]}";
do
	echo "This is what will be downloaded to make the mash database: $val"
  	valUp="${val:1:-1}"				## Remove quotes
  	valUp=${val//[[:blank:]]/}		## Remove space

	echo "Beginning to dowload genomes from NCBI..."

		if [[ -z "$assembly" ]] ; then
			echo "Assembly level is not specified..."

		## With datasets version change 12.2.0 -> 15.2.0
		## complete_genome -> complete & no longer required to specify files to exclude
			datasets download genome taxon "$val" --dehydrated --assembly-source genbank \
--filename $valUp.zip --assembly-level complete,chromosome,scaffold,contig 

		elif [[ -n "$assembly"  ]]; then
			echo "Assembly level specified. Downloading only $assembly ..."
			datasets download genome taxon "$val" --dehydrated --assembly-source genbank \
 --filename $valUp.zip --assembly-level "$assembly"   2>/dev/null || error_handler_assembly
		fi

		## Confirm zip file avaiable
		count=$(ls -1 *.zip 2>/dev/null | wc -l)
		if [[ $count -ne 0 ]]; then
			echo "Checking for zip files, should be present after downloading data..."
			echo "Number of zip files: $count"
			echo "Unzipping the associated files..."
			unzip $valUp.zip -d $valUp

			datasets rehydrate --directory $valUp

			#shorten path 
			speciesdownload=$valUp"/ncbi_dataset/data"

			# Check for genomes with NCBI taxonomy issues and add GenBank/RefSeq IDs to the exclusion list
			# If a genome GCA ID is added to the list, it will be deleted 
			awk '{if (!/OK/) print $1}' "$subfolder/$speciesdownload/assembly_data_report.jsonl" | \
				grep -o "GCA_..........." >> excluded_genomes.tmp

			## Get 'unculture' legionella species GCA ids and adde to a list (excluded_genomes)
			awk '{if(/uncultured/) print $1}' $subfolder/$speciesdownload/assembly_data_report.jsonl | \
				grep -o "GCA_..........." >> excluded_genomes.tmp

			## Exclude genomes with completeness below 93.00 (range 0 - 100). 
			## Inconsistencies exist in the assembly file fields between isolates; for example,
			## some lack contamination estimates. 
			## This script retrieves data from genomes with completeness ≥ 93.00 
			## (printing fields $1 and $3) and those without (printing fields $1 and $4).
			## Genome GCA IDs that are excluded are added to the excluded_genomes list.
			
			cat $subfolder/$speciesdownload/assembly_data_report.jsonl |  grep -o "completeness.*" | grep -o ".*organism" | \
				awk -F , '{ print $1 $3 }' | grep -v "contamination" | awk -F\" '{ print $2 " " $5 }' | awk -F : '{ print $2 }' | \
				awk '{ if( $1 < 93.00) print $2 }' >> excluded_genomes.tmp
	
			cat $subfolder/$speciesdownload/assembly_data_report.jsonl | grep -o "completeness.*" | grep -o ".*organism" | \
				awk -F , '{ print $1 $4 }' | grep -v "contamination" | awk -F\" '{ print $2 " " $5 }' | awk -F : '{ print $2 }' | \
				awk '{ if( $1 < 93.00) print $1 " " $2 }' | grep GCA | awk '{ print $2 }' >> excluded_genomes.tmp

			cat excluded_genomes.tmp | uniq -u >> excluded_genomes.txt

			## Now remove folders for genomes listed in exclusion file
			TO_BE_DEL="excluded_genomes.txt"
			while read -r file ; do
 
				rm -r $speciesdownload/"$file" > /dev/null 2>> error.log
			done < "$TO_BE_DEL"

			cp excluded_genomes.txt $basefolder
			rm excluded_genomes.tmp

			cd $valUp/ncbi_dataset/data

			## Check for plasmids and remove
			## Get fasta header /^>/ & assign header with plasmid to p
			## Remove plasmid from genome by specifying not P (!p)
			## When it is in the same file as FASTA data, which happens
			## data is concatenated together
			
			echo "Checking for .fna files..."
			files=( */*.fna )

			if [ ${#files[@]} -eq 0 ]; then
				echo "No .fna files found. Creating a placeholder file." >> temp_file-noFNA.fna
			else
				for i in "${files[@]}"; do
					if [ -f "$i" ]; then  # Check if it's a file
						echo "Checking for plasmids in $i. This can take some time..."
						awk '/^>/ { p = ($0 ~ /plasmid/) } !p' "$i" > "${i%.*}_cleaned.fna"
						mv "${i%.*}_cleaned.fna" "$i"
					else
						echo "$i is not a file, skipping..." > /dev/null 2>> error.log
					fi
				done
			fi

			# For each file that matches the criteria (-name and if a file is just a plasmid), then remove
			find . \( -name "chrunnamed*.unlocalized.scaf.fna" -o -name "cds_from_genomic.fna" -o -size 0 -type f \) -exec rm -rf {} +
			find . -name "*.fna" -exec grep -l "plasmid" {} \; -exec rm {} +

			## Make summary file of the downloaded data
			## Get opposite of grep, so do not pass excluded genome information for file manipulation/checking
			echo "Making $val map file for file name conversion and converting file names..."

			dataformat tsv genome --inputfile *.jsonl --fields organism-name,accession,assminfo-paired-assmaccession --force  | \
			grep -vFwf $subfolder/excluded_genomes.txt  | \
			awk 'FNR==1 { header = $0; print }  $0 != header' | \
			sed 's/\//-/g' | \
			sed 's/ /_/g' | tee temp | \
			awk '{ print $1 "_" $2 }' > temp2 				## Now combine column 1 with underscore and column 2

			## Remove headers
			sed -i '1d' temp 
			sed -i '1d' temp2 

			## Make final map file
			cut -f2 temp | paste -d " " temp2 - | \
			awk '{ print $2 ".fna" " " $1 }' > mapFinal$valUp.txt 				## File w/2 columns: old file name & new file name

			## Rename *.fna to folderName.fna to handle unplaced scaffolds and duplicate file names between isolates
			for d in */
			do
				filepath=$d*.fna 
				mv $filepath "$(dirname "$filepath")/$(dirname "$filepath").fna" > /dev/null 2>> error.log
			done
	
			for file in */*.fna; do
			## Check if the file does not end with -noFNA.fna
			## These files are created, if after exclusion there is no .fna files to continue process
				if [[ ! "$file" == *"-noFNA.fna" ]]; then
					cp "$file" "$subfolder/allDownload" > /dev/null 2>> error.log
				else 
					cp "$file" "$basefolder" 	
				fi
			done
			
			## Rename files using mapfile
			cd $subfolder/allDownload 
			awk -F " " 'system("mv " $1 " " $2 ".fna")'  $subfolder/$valUp/ncbi_dataset/data/mapFinal$valUp.txt
			
			## Move back to base directory of genomesDownloaded_timestamp
			cd $subfolder
			echo "Summarizing the entire download..."

####################
##CREATE SUMMARIES##
####################
## INPUT TSV file form NCBI command line tool dataformat 
## OUTPUT TSV file with species and genebank accession downloaded
		
			dataformat tsv genome --package $valUp.zip --fields organism-name,accession,assminfo-paired-assmaccession --force | \
			grep -vFwf $subfolder/excluded_genomes.txt | \
			awk 'FNR==1 { header = $0; print }  $0 != header' | \
			grep -v 'Legionella sp\.' | \
			grep -v 'Legionella sp\. ' >> downloadedData.tsv    ## Must include space before or get rid of Lp unknown species
	
			## Create a txt file of a count of all species downloaded
			cat downloadedData.tsv | sed 1d | cut -f1 | cut -f2 -d ' ' | sort | uniq -c > speciesCount.txt
	
			## Get total number of isolates; store and compare to those in allDownload directory
			countFile=$(awk '{ SUM+=$1}END{print SUM }' speciesCount.txt)
			awk '{ SUM+=$1}END{print SUM " Total Isolates"}' speciesCount.txt >> speciesCount.txt
	
			cd allDownload || { echo "Failed to change directory to allDownload"; exit 1; }
	
			## Exclude legionella that is not identified to species and endosymbionts
			if [[ "${species^}" == "Legionella" ]]; then
				echo "Removing Legionella only identified to genus (e.g., Legionella sp.)..."
				rm Legionella_sp._*.fna
			else
				echo "No extra files to remove for non-Legionella species."
			fi
	
			## Count number of files in allDownload folder with those in speicesCount file for comparison
			countFolder=$(ls | wc -l)

			if [[ $countFile  -eq  $countFolder ]]; then
				echo "The number of isolates in the speciesCount file matches the number of files in allDownload directory.";
			else
				echo "Mismatch between speciesCount file and downloaded files. Investigate or retry the script.";
				exit 1
			fi
	
			## Move files up to basefolder to allow easier copying via nextflow process
			count=$(ls -1 *.fna 2>/dev/null | wc -l)
			
			if [[ $count -gt 0 ]]; then
				mv *.fna "$basefolder" || echo "Failed to move .fna files. Exiting"    
			else
				echo "No .fna files generated. This is a placeholder. " >> $basefolder/$valUp-noFNA.fna
			fi
		else
			echo "Exiting the program."
			echo "----------------------------------------"
		fi
done