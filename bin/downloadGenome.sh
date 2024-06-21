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
         echo 'Invalid Option. Exiting.' >&2 
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
      if [ $count != 0 ]; then

        echo $count
        echo "Unzipping the associated files..."
        unzip $valUp.zip -d $valUp
        echo " "

        datasets rehydrate --directory $valUp

        ## Checking for genomes where NCBI taxonomy is NOT OK, adding the genebank/refseq id to a list, 
        ## if found and then deleting those folders via while loop that reads the created exclude_genomes file
        cat $basefolder/genomesDownloaded_$timestamp/$valUp/ncbi_dataset/data/assembly_data_report.jsonl | awk '{if (!/OK/) print $1}' | grep -o "GCA_..........." >> excluded_genomes.tmp
        
        ## remove genomes without a completeness estimate
        #cat $basefolder/genomesDownloaded_$timestamp/$valUp/ncbi_dataset/data/assembly_data_report.jsonl | awk '{ if(!/"completeness"/) print $1 }' | grep -o "GCA_......." 
        
        ## get 'unculture' legionella species and remove those genomes
        ## If files are already listed to be removed in a previous command than an error message gets printed, but this is not an error, it is just that it can not find the already deteled files.
        cat $basefolder/genomesDownloaded_$timestamp/$valUp/ncbi_dataset/data/assembly_data_report.jsonl |  awk '{if(/uncultured/) print $1}' | grep -o "GCA_..........." >> excluded_genomes.tmp
        
        ## I am sure there is a better way to do this, but have not found it. I want to exclude genomes which 
        ## have a lower level of completeness (< 93.00), but all fields in the assembly file are not consistent, 
        ## meaning that some do not have a contamiation estimate, so I get the data which does and then parse 
        ## it here (print $1 $3) and below I get the data that does not and parse it (print $1 $4)
        cat $basefolder/genomesDownloaded_$timestamp/$valUp/ncbi_dataset/data/assembly_data_report.jsonl |  grep -o "completeness.*" | grep -o ".*organism" | awk -F , '{ print $1 $3 }' | grep -v "contamination" | awk -F\" '{ print $2 " " $5 }' | awk -F : '{ print $2 }' | awk '{ if( $1 < 93.00) print $2 }' >> excluded_genomes.tmp
 
        cat $basefolder/genomesDownloaded_$timestamp/$valUp/ncbi_dataset/data/assembly_data_report.jsonl | grep -o "completeness.*" | grep -o ".*organism" | awk -F , '{ print $1 $4 }' | grep -v "contamination" | awk -F\" '{ print $2 " " $5 }' | awk -F : '{ print $2 }' | awk '{ if( $1 < 93.00) print $1 " " $2 }' | grep GCA | awk '{ print $2 }' >> excluded_genomes.tmp

        cat excluded_genomes.tmp | uniq -u >> excluded_genomes.txt
 
        TO_BE_DEL="excluded_genomes.txt"
        while read -r file ; do
 
          rm -r $valUp/ncbi_dataset/data/"$file" >/dev/null 
        done < "$TO_BE_DEL"

        cp excluded_genomes.txt $basefolder
        rm excluded_genomes.tmp
        
        cd $valUp/ncbi_dataset/data
          
        ## Check for plasmids and remove
        echo "Checking for plasmids. This can take a some time..."
        for i in */*.fna
        do
          awk '/^>/ { p = ($0 ~ /plasmid/) } !p' $i > ${i%\.*}_cleaned.fna
          mv ${i%\.*}_cleaned.fna $i
        done

        find . -name "chrunnamed*.unlocalized.scaf.fna" -exec rm -rf {} \;        ## These are plasmid files also
        find . -name "*.fna" -exec grep "plasmid" {} \; -exec rm {} \;
        find . -name "cds_from_genomic.fna" -exec rm -rf {} \;                    ## Remove these files. Downloaded via conda
        find . -size 0 -type f -delete                                            ## Remove files with zero bytes

        ## Make summary file of the downloaded data
        echo "Making $valUp map file for file name conversion and converting file names..."

        dataformat tsv genome --inputfile *.jsonl --fields organism-name,accession,assminfo-paired-assmaccession --force >> temp
        
        ## Get opposite of grep, so do not pass excluded genome information for file manipulation/checking because deleted
        grep -vFwf $basefolder/genomesDownloaded_$timestamp/excluded_genomes.txt temp > temp2
        
        awk 'FNR==1 { header = $0; print }  $0 != header' temp2 > temp3 #downloaded-$valUp.tsv ## Remove duplicate header
        
        sed -i 's/\//-/g' temp3  ##downloaded-$valUp.tsv

        ## Replaces spaces with dash
        sed -i 's/ /_/g' temp3  ##downloaded-$valUp.tsv

        ## Now combine column 1 with underscore and column 2
        awk '{ print $1 "_" $2 }' temp3 > temp4 #   downloaded-$valUp.tsv > map2$valUp.txt

        ## Remove headers
        sed -i '1d' temp3 ##downloaded-$valUp.tsv
        sed -i '1d' temp4 ##map2$valUp.txt
        #cp temp3  downloaded-$valUp.txt

        ## Make final map file
        cut -f2 temp3 | paste -d " " temp4 - > temp5

        ## Change *.fna to only folderName.fna. This deals with unplaced scaffolds and
        ## File names that are duplicated between isolates
        for d in */
        do
          FILEPATH=$d*.fna 
          mv $FILEPATH "$(dirname "$FILEPATH")/$(dirname "$FILEPATH").fna"
        done

        ## Makes file of two columns old file name and new file name
        awk '{ print $2 ".fna" " " $1 }' temp5 > mapFinal$valUp.txt

        ## Move to common folder
        mkdir common
        
        cp */*.fna common
        
        cp mapFinal$valUp.txt common

        ## Rename files using mapfile
        cd common
        awk -F " " 'system("mv " $1 " " $2 ".fna")'  mapFinal$valUp.txt

        ## Move all converted *.fna files from species common to alldownload
        mv *.fna $basefolder/genomesDownloaded_$timestamp/allDownload
        cd ..
       # rm -r common 
       # rm temp temp2 temp3 temp4 temp5 mapFinal$valUp.txt
        
        ## Move back to base directory of genomesDownloaded_timestamp
        cd $subfolder
        echo "Summarizing the entire download..."

####################
##CREATE SUMMARIES##
####################
## This uses NCBI command line tool dataformat on the zipped files
## Output is a tsv file with species and genebank accession downloaded

        
        dataformat tsv genome --package $valUp.zip --fields organism-name,accession,assminfo-paired-assmaccession --force >> temp1
        grep -vFwf $basefolder/genomesDownloaded_$timestamp/excluded_genomes.txt temp1 >> temp2

        awk 'FNR==1 { header = $0; print }  $0 != header' temp2 > temp3    ## Remove duplicate header if doing multiple species

        grep -v 'Legionella sp\. ' temp3 > temp4
        grep -v 'Legionella sp\.' temp4 > downloadedData.tsv    ## Must include space before or get rid of Lp unknown species

        ## Create a txt file of a count of all species downloaded
        cat downloadedData.tsv | sed 1d | cut -f1 | cut -f2 -d ' ' | sort | uniq -c > speciesCount.txt

        ## Get total number of isolates; store and compare to those in allDownload directory
        countFile=$(awk '{ SUM+=$1}END{print SUM }' speciesCount.txt)

        awk '{ SUM+=$1}END{print SUM " Total Isolates"}' speciesCount.txt >> speciesCount.txt

        ## Remove temp file within base directory genomesDownloaded_timestamp
        rm temp1 temp2 temp3 temp4 
        

#################################
##SECOND FILE CLEANUP AND CHECK##
#################################
        cd allDownload

        ## Exclude legionella that is not identified to species and endosymbionts
        if [[ "${species^}" == "Legionella" ]]; then
          echo "Removing Legionella files where the isolate is not identified to a recognized species..."
          rm Legionella_sp._*.fna
        else
          echo "This was not Legionella identified to only to sp. (Legionella sp.). Thus, no extra files to remove..."
        fi

        ## Count number of files in folder with those in speicesCount file for comparison
        countFolder=$(ls | wc -l)

        if [[ $countFile  -eq  $countFolder ]]; then
          echo "Good, the number of isolates in the speciesCount file matches the \
number of files that were copied to the allDownload directory.";
        else
          echo "Hmm, the number of isolates in the speciesCount file does not match \
the number of files that were copied to the allDownload directory. You \
#can either investigate the output/error files or just try running the script\
#again as sometimes there are communication issues between HPC and NCBI.";
        exit 1
        fi

        ## Move files up to basefolder to all easier copying via nextflow process
        count=`ls -1 *.fna 2>/dev/null | wc -l`
        if [ $count != 0 ]; then
          mv *.fna $basefolder ## should be mv     
          rm -rf $basefolder/genomesDownloaded_$timestamp/$valUp
        else
           echo "There were no .fna files generated"
        fi
        echo "Exiting the program."
        echo "-----------------------------"
  fi
done
