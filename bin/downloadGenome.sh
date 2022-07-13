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
   echo "Download all genomes from NCBI using their command line tools (datasets/dataformat). \
Then convert GCA#.fna files to be named Genus_species_GCA#.fna"
   echo
   echo "Syntax: downloadGenome.sh [-c|s|h]"
   echo "Options:"
   echo " -c    Required: Activate the conda enviroment? (T/F). \
Assumes conda environment named ncbi_datasets"
   echo " -s    Required: Genus or species to download. Can use multiple -s flags. \
Example: -s \"legionella\" -s \"tatlockia micdadei\" "
   echo " -h    Print this help."
   echo " "
   echo "Example command: downloadGenome.sh -c F -s \"legionella\" "
   echo " "
}

## Get the options
while getopts ":c:s:h" option; do
   case $option in
      h) ## Display Help
         Help
         exit;;
      c) ## Activate the conda env? (T/F)
         conda=$OPTARG;;
      s) ## Specifiy the species you want to download. Can be multiple, seperated by a space. If using genus species then place in quotes ("tatlockia micdadei")
         species+=("$OPTARG");;
   esac
done

###################
##PARSE THE FLAGS##
###################

## Make sure both -c and -s flags are included and not empty
if [[ $conda == "" && $species == "" ]]
then
    echo 'You must supply both the -c and -s flags along with their required request.
Please look at the help menu. Exiting.'
    echo ""
    exit 1
elif [[ $conda == "" ]]
then
    echo 'You did not provide any information after the -c flag (T/F). Exiting.'
    exit 1
elif [[ $species == "" ]]
then
    echo 'You did not provide an organism to download for the -s flag. Exiting.'
    exit 1
fi

## Determine what to do based on input from -c (conda) flag
if [[ $conda == @(True|true|T|t) && $species ]]
then
    echo 'Activating the conda environment, assuming it is called ncbi_datasets...'
    eval "$(conda shell.bash hook)"
    conda activate ncbi_datasets
    condaAct=`echo $CONDA_DEFAULT_ENV`
    echo "This is your conda enviroment that is activated:" $condaAct
    if [[ $condaAct == 'ncbi_datasets' ]]
    then
        echo 'These are the tools in your conda enviroment....'
        conda list -n ncbi_datasets
    else
        echo "This tool assumes the ncbi datasets cli conda environment is called ncbi_datasets. Exiting."
        exit 1
    fi
elif [[ $conda == @(False|false|F|f) && $species ]]
then
    echo "Confirming both NCBI datasets and dataformat tools are available..."
  ## Check that both tools are available. If not then exit
    command -v dataformat >/dev/null 2>&1 || { echo >&2 "NCBI dataformat is not installed.  Exiting."; exit 1; }
    command -v datasets >/dev/null 2>&1 || { echo >&2 "NCBI datasets is not installed.  Exiting."; exit 1; }
    #probably should have a response that the tools are avialble
    echo "Great both tools available to access NCBI..."
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
    echo "Good, a $FILE summary file does not already exist, will continue..."
fi

FILE2=speciesCount.txt
if [[ -f "$FILE2" ]]; then
    echo "$FILE2 exists, this script will overwrite to previous $FILE2 \
and you will lose that summary file. Exiting"
    exit 1
else
    echo "Good, the $FILE2 summary file doesn't already exist, will continue..."
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

for val in "${species[@]}";
do

  valUp="${val:1:-1}"                                                           ## Remove quotes
  valUp=${val//[[:blank:]]/}                                                    ## Remove space

  echo "Beginning to dowload genomes from NCBI..."

  datasets download genome taxon "$val" --dehydrated --exclude-gff3 \
  --exclude-protein --exclude-rna --assembly-source genbank \
  --filename $valUp.zip

  if [[ $? -ne 0 ]] ; then
      echo "It appears that the name you specified via -s flag is not recognized by NCBI, please check spelling. Exiting."
      exit 1
  else
      echo "Good, the names are recognized by NCBI. Continuing..."
  fi

#TO DO: add better documentation here:
 #added this if statement to deal with of unzipping.
  if test -z "$CONDA_DEFAULT_ENV"; then
     echo "Okay this is when -c F:" $condaAct
     unzip $valUp.zip -d $valUp
  else
     echo "Okay this is when -c T:" $condaAct
     unzip $valUp.zip -d $valUp
     #7z x $valUp.zip -o*

  fi

  datasets rehydrate --directory $valUp

  cd $valUp/ncbi_dataset/data

  ## Check for plasmids and remove
  echo "Checking for plasmids..."

  for i in */*.fna
  do
    awk '/^>/ { p = ($0 ~ /plasmid/) } !p' $i > ${i%\.*}_cleaned.fna
    mv ${i%\.*}_cleaned.fna $i
  done

  find . -name "chrunnamed*.unlocalized.scaf.fna" -exec rm -rf {} \;            ## These are plasmid files also
  find . -name "*.fna" -exec grep "plasmid" {} \; -exec rm {} \;
  find . -name "cds_from_genomic.fna" -exec rm -rf {} \;                        ## Remove these files. Downloaded via conda
  find . -size 0 -type f -delete                                                ## Remove files with zero bytes

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

done

echo "Summarizing the entire download..."

####################
##CREATE SUMMARIES##
####################
## This uses NCBI command line tool dataformat On the zipped files
## Output is a tsv file with species and genebank accession downloaded

for val in "${species[@]}"
do
  valUp="${val:1:-1}"                                                           ## Remove quotes
  valUp=${val//[[:blank:]]/}                                                    ## Remove space

  dataformat tsv genome --package $valUp.zip \
  --fields organism-name,assminfo-genbank-assm-accession,assminfo-refseq-assm-accession >> temp
done

awk 'FNR==1 { header = $0; print }  $0 != header' temp > downloadedData.tsv     ## Remove duplicate header if doing multiple species

## Exclude legionella that is not identified to species or is endosymbionts.
## If other species, will give error of no such file. Could ignore but better
## To not unless for legionella....

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
  echo "Removing Legionella endosymbiont files and files where the isolate is only identified to genus..."
  rm Legionella_sp._*.fna
  rm Legionella_endosymbiont*.fna
else
  echo "This was not either Legionella endosymbiont or those identified to species (Legionella sp.). Thus, no extra files to remove..."
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
rm -rf $basefolder/genomesDownloaded_$timestamp/allDownload
rm -rf $basefolder/genomesDownloaded_$timestamp/$valUp

echo "Exiting the program."
echo " "
