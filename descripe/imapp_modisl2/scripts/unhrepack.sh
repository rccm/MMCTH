#!/bin/bash

# Script to uncompress a HDF file which may have been
# internally compressed with hrepack (e.g., MODIS L1B files from LADS)
# Translated from csh to bash by R. Cintineo July 2013

echo
echo
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Uncompressing files"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo

# Check arguments
if [ $# == 0 ]; then
  echo "Usage: unhrepack.sh file.hdf"
  echo "or"
  echo "       unhrepack.sh *.hdf"
  exit 1
fi

# Loop through each file in the argument list
while (( $# ))
do

  # Set input and output file names
  input_file=$1
  output_file=$1.unhrepack
  
  # Uncompress this file
  hrepack -i $input_file -o $output_file -t '*:NONE'
  if [ $? != 0 ]; then
    echo "Could not uncompress input file "$input_file
    /bin/rm $output_file
    shift
    continue 
  fi
  
  # Replace input file with uncompressed output file
  /bin/mv $output_file $input_file
  if [ $? != 0 ]; then
    echo "Could not replace input file "$input_file
    rm $output_file
  else
    echo "Uncompressed "$input_file
  fi
  
  # Get next argument
  shift
  
done

exit 0
