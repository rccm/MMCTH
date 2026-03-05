#!/bin/bash

# Script to internally compress HDF file(s) with hrepack
# Translated from csh to bash by R. Cintineo July 2013

# Check arguments
if [ $# == 0 ]; then
  echo "Usage: hrepack.sh file.hdf"
  echo "or"
  echo "       hrepack.sh *.hdf"
  exit 1
fi

# Loop through each file in the argument list
while (( $# ))
do

  # Set input and output file names
  input_file=$1
  output_file=$1.hrepack
  
  # Compress this file
  hrepack -i $input_file -o $output_file -t '*:GZIP 1'
  if [ $? != 0 ]; then
    echo "Could not compress input file "$input_file
    /bin/rm $output_file
    shift 
  fi
  
  # Replace input file with compressed output file
  /bin/mv $output_file $input_file
  if [ $? != 0 ]; then
    echo "Could not replace input file "$input_file
    rm $output_file
  else
    echo "Compressed "$input_file
  fi
  
  # Get next argument
  shift
  
done
exit 0
