#!/bin//bash

files=$(ls | wc -l)

echo "Running mash sketch on $files files"
for file in *.fna
do
  mash sketch $file -k 25 -s 100000 -p 8
done

timestamp=$(date +%d-%m-%Y)
echo "Completed generating individual mash files. Now will generate Mash database."
mash sketch -o MASH_Legionella_master_sketch_$timestamp.msh *.msh

echo 'Completed making the mash database. Exiting'
