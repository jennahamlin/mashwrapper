#!/bin/bash

## Author J. Hamlin
## Date 01/05/2022
## This script is to download all available *legionella* and
## *legionella*-associated genomes from NCBI using the NCBI datasets command
## line tools. The NCBI dataformat command line tool then summarizes the output
## by species count. Addiitonally, it converts non-species descript file names
## (e.g., GCA_XXXX) to include genus and species in the file name. In theory,
## this tool should download data for any species located on NCBI by specifying
## the organism of interest through the -s flag. Though I have only tested the
## tool for the species listed below so careful double checking of download and
## file conversion should be done when used for other species.

## species tested on
## "legionella" "fluoribacter" "tatlockia micdadei"

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
   echo "Syntax: downloadGenome.sh [-c|s|h]"
   echo "Options:"
   echo " -a   Optional: Assembly-level Restrict assemblies to be \
chromosome, complete_genome, contig, scaffold"
   echo " -c   Required: Activate the conda enviroment? (T/F). \
Assumes conda environment named ncbi_datasets"
   echo " -s   Required: Download genus or species. Can use multiple -s flags. \
Example: -s \"legionella\" -s \"tatlockia micdadei\" "
   echo " -h   Print this help."
   echo " "
   echo "Example command: downloadGenome.sh -c F -s \"legionella pneumophila\" -a complete_genome "
   echo " "
}

## Get the options
while getopts ":a:c:s:h" option; do
   case $option in
      h) ## Display Help
         Help
         exit;;
      a) ## Assembly-level
         assembly=$OPTARG;;
      c) ## Activate the conda env? (T/F)
         conda=$OPTARG;;
      s) ## Specifiy the organism to download. Can be multiple, seperated by a
         ## space. Using genus species, place in quotes ("tatlockia micdadei")
         species+=("$OPTARG");;
   esac
done

###################
##PARSE THE FLAGS##
###################

## Make sure both -c and -s flags are included and not empty
if [[ $conda == "" && $species == "" ]]
then
    echo 'Please supply both -c and -s flags along with their required request.
Please look at the help menu, if you need additional information. Exiting.'
    echo ""
    exit 1
elif [[ $conda == "" ]]
then
    echo 'You did not provide any information after the -c flag (T/F).
Please look at the help menu, if you need additional information. Exiting.'
    exit 1
elif [[ $species == "" ]]
then
    echo 'You did not provide an organism to download for the -s flag. Exiting.
Please look at the help menu, if you need additional information. '
    exit 1
fi

## Determine what to do based on input from -c (conda) flag
## First check if conda is specifed as True
if [[ $conda == @(True|true|T|t) && $species ]]
then
  ## Determine if conda == T request is from nextflow or not, by that I mean
  ## is the conda environment loaded  locally or requested to be used via
  ## nextflow. $nf exported from the downloadGeome.nf module, if from nextflow.
    if [[ -z "$nf" ]]
    then
        echo $nf
        echo 'The variable check to determine if next flow is running is empty.
Activating your local conda environment. Continuing...'
        eval "$(conda shell.bash hook)"
        conda activate ncbi_datasets
        condaAct=`echo $CONDA_DEFAULT_ENV`
        echo "This is your local conda enviroment that is activated:" $condaAct
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
  echo "Confirming both NCBI datasets/dataformat tools and Mash are available..."
      ## Check that both tools are available. If not then exit
  command -v dataformat >/dev/null 2>&1 || { echo >&2 "NCBI dataformat is not installed.  Exiting."; exit 1; }
  command -v datasets >/dev/null 2>&1 || { echo >&2 "NCBI datasets is not installed.  Exiting."; exit 1; }
  #command -v mash >/dev/null 2>&1 || { echo >&2 "Mash is not installed.  Exiting."; exit 1; }
      ## Inform user that the tools are accessible for when conda is false
  echo "Great tools available to access NCBI and run Mash. Continuing.."
else
  echo 'Unable to activate a conda environment, find the tools in a bin folder,
or confirm that the tool is using a container. Exiting.'
  exit 1
fi

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

####################
##DOWNLOAD GENOMES##
####################
## Loop through species array and download genomes
## Can change source to genbank (GCA) or refseq (GCF)
## I use genebank option as this has a larger number of genomes

# echo "If assembly-level was specified, then this was the level of restriction \
# for genomes to download: $assembly"

error_handler()
{
  #echo 'This will exit from this isolate but test the others...'
  echo "No $assembly" files available. Creating a file place holder for this species: $valUp. Exiting.
  cd ..
  echo "There are no $assembly files avilable at the level specified. Exiting." > $valUp-$assembly-noFNA.fna
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
        datasets download genome taxon "$val" --dehydrated --exclude-gff3 \
--exclude-protein --exclude-rna --assembly-source genbank \
--filename $valUp.zip --assembly-level complete_genome,chromosome,scaffold,contig

      elif [[ -n "$assembly"  ]]; then
        echo "Assembly level is specified and will only download $assembly..."
        datasets download genome taxon "$val" --dehydrated --exclude-gff3 \
 --exclude-protein --exclude-rna --assembly-source genbank \
 --filename $valUp.zip --assembly-level "$assembly"   2>/dev/null || error_handler
      fi

## When running testGet or get_database with singularity, then unzip will complain. Error = Bad length
## But process still runs to completion and is successful. As far as I can tell, it maybe an issue
## with the unzip command which is in the BusyBox instance associated with singularity. When looking
## at the unzip github (https://github.com/brgl/busybox/blob/master/archival/unzip.c), a comment is
## included which says "Don't die. Who knows maybe len calculation was botched somewhere. After all,
## crc matched!" This is associated with the validate comparsion if statement (Line 392). Using
## `unzip -l` within the continue provides an estimate of length, but not in human readable form and
## does not match if ls -lh after unzipping.

## TO DO: check this on another hpc
##if [[  -z "$CONDA_DEFAULT_ENV" ]]; then
##   conda default env should be empty if using a container
##   echo "This is when -c is False as in when a container is used. No conda environment should be listed:" $condaAct
      count=`ls -1 *.zip 2>/dev/null | wc -l`
  #if [[ -f $valUp.zip ]]; then
      if [ $count != 0 ]; then
        #echo "Line 221 - and this is count"
        #echo $count
        echo "Unzipping the associated files..."
        unzip $valUp.zip -d $valUp
        echo " "
        echo 'NOTE TO USER: unzip: bad length is nothing to worry about. Tool runs
to completion successfully. It might be with a len calculation with 
unzip in the BusyBox instance associated with the container.
See: https://github.com/brgl/busybox/blob/master/archival/unzip.c '
        echo " "
#elif [[ $condaAct != "ncbi_datasets" ]]; then
#     echo "This is when -c is False w/o config file. Assume running w/modules loaded. No conda env should be list:" $condaAct
#     unzip $valUp.zip -d $valUp
#elif ! [[ -x "$(command -v 7z)" ]]; then
#     echo "Checking for the program 7z..."
#     echo "7z installs with conda loaded from nf module"
#     unzip $valUp.zip -d $valUp
#else
#     echo "This is when -c is True as the user specified using the conda environment:" $condaAct
#     7z x $valup.zip -o*
#fi
        datasets rehydrate --directory $valUp

        cd $valUp/ncbi_dataset/data

  ## Check for plasmids and remove
        echo "Checking for plasmids..."
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
        echo "Making $valUp map file for file name conversion..."

        dataformat tsv genome --inputfile *.jsonl \
--fields organism-name,assminfo-genbank-assm-accession,assminfo-refseq-assm-accession >> temp

        awk 'FNR==1 { header = $0; print }  $0 != header' temp > downloaded-$valUp.tsv ## Remove duplicate header
        rm temp

        sed -i 's/\//-/g' downloaded-$valUp.tsv

        ## Replaces spaces with dash
        sed 's/ /_/g' downloaded-$valUp.tsv > map1$valUp.txt

        ## Now combine column 1 with underscore and column 2
        awk '{ print $1 "_" $2 }' map1$valUp.txt > map2$valUp.txt

        ## Remove headers
        sed -i '1d' map1$valUp.txt
        sed -i '1d' map2$valUp.txt

        ## Make final map file
        cut -f2 map1$valUp.txt | paste -d " " map2$valUp.txt - > mapFinal$valUp.txt

        ## Change *.fna to only folderName.fna. This deals with unplaced scaffolds and
        ## File names that are duplicated between isolates
        for d in */
        do
          FILEPATH=$d*.fna
          mv $FILEPATH "$(dirname "$FILEPATH")/$(dirname "$FILEPATH").fna"
        done

        ## Makes file of two columns old file name and new file name
        awk '{ print $2 ".fna" " " $1}' mapFinal$valUp.txt > mapFinal2$valUp.txt

        ## Move to common folder
        mkdir common
        cp */*.fna common
        cp mapFinal2$valUp.txt common

        ## Rename files using mapfile
        cd common
        awk -F " " 'system("mv " $1 " " $2 ".fna")'  mapFinal2$valUp.txt

        ## Move all converted *.fna files from species common to alldownload
        cp *.fna $basefolder/genomesDownloaded_$timestamp/allDownload

        ## Move out of common folder
        cd ..
        rm -r common

        ## Move back to base directory of genomesDownloaded_timestamp
        cd $subfolder
        echo "Summarizing the entire download..."

####################
##CREATE SUMMARIES##
####################
## This uses NCBI command line tool dataformat on the zipped files
## Output is a tsv file with species and genebank accession downloaded

        dataformat tsv genome --package $valUp.zip \
--fields organism-name,assminfo-genbank-assm-accession,assminfo-refseq-assm-accession >> temp
#done

        awk 'FNR==1 { header = $0; print }  $0 != header' temp > downloadedData.tsv     ## Remove duplicate header if doing multiple species

        ## Exclude legionella that is not identified to species or is endosymbionts.

        grep -v 'Legionella endosymbiont' downloadedData.tsv > temp
        grep -v 'Legionella sp. ' temp > downloadedData.tsv

        ## Create a txt file of a count of all species downloaded
        cat downloadedData.tsv | sed 1d | cut -f1 | cut -f2 -d ' ' | sort |uniq -c > speciesCount.txt

        ## Get total number of isolates; store and compare to those in allDownload directory
        countFile=$(awk '{SUM+=$1}END{print SUM }' speciesCount.txt)

        awk '{SUM+=$1}END{print SUM " Total Isolates"}' speciesCount.txt >> speciesCount.txt

        ## Remove temp file within base directory genomesDownloaded_timestamp
        rm temp

#################################
##SECOND FILE CLEANUP AND CHECK##
#################################
        cd allDownload

## Exclude legionella that is not identified to species and endosymbionts
        if [[ "$species" == "legionella" ]]; then
          echo "Removing Legionella endosymbiont files and files where the isolate is
only identified to genus..."
          rm Legionella_sp._*.fna
          rm Legionella_endosymbiont*.fna
          rm uncultered_Legionella*.fna
          rm Legionella genomosp.*
        else
          echo "This was not either Legionella endosymbiont or those identified to
species (Legionella sp.). Thus, no extra files to remove..."
        fi

## Count number of files in folder with those in speicesCount file for comparison
        countFolder=$(ls | wc -l)

        if [[ $countFile  -eq  $countFolder ]]; then
          echo "Good, the number of isolates in the speciesCount file matches the \
number of files that were copied to the allDownload directory.";
        else
          echo "Hmm, the number of isolates in the speciesCount file does not match \
the number of files that were copied to the allDownload directory. You \
can either investigate the downloadG.e***/downloadG.o*** file or just try \
running the script again as sometimes there are communication issues \
between HPC and NCBI.";
        fi
## Move files up to basefolder to all easier copying via nextflow process
        mv *.fna $basefolder
       # rm -rf $basefolder/genomesDownloaded_$timestamp/allDownload
        rm -rf $basefolder/genomesDownloaded_$timestamp/$valUp

        echo "Exiting the program."
        echo "-----------------------------"
  fi
done
