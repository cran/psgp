#!/bin/bash
# This script should be run once to set up the folder structure.
# It merely creates symbolic links to the relevant folders in 
# gptk/libgptk.

# Folder where the libgptk source is installed
R_SRC=`pwd`
GPTK_ROOT=`pwd | sed -e "s/\/R\/psgp\/src//"`
GPTK_SRC=$GPTK_ROOT/libgptk

# FILES/FOLDERS TO LINK TO
LINK_TO=("Makevars.in"
         "astonGeostats.cpp"
         "covarianceFunctions"
         "io"
         "likelihoodModels"
	 "RInterface.cpp"
	 "astonGeostats.h"
	 "gaussianProcesses"
	 "itppext"
	 "optimisation")

# LINK NAMES
LINK_NAME=( "${LINK_TO[@]}" )

# NUMBER OF LINKS
FILE_COUNT=${#LINK_TO[*]}

echo "\n** This script will set up the correct structure for"
echo "** the R package. It only needs be run once."

echo "- PSGP R package source folder: " $R_SRC
echo "- libgptk source folder: " $GPTK_SRC

# DELETE EXISTING LINKS
echo
echo "** Deleting existing symbolic links"
for i in `seq 0 $(( FILE_COUNT-1 ))`
do
  echo Deleting ${LINK_NAME[$i]}
  rm ${LINK_NAME[$i]}
done

# CREATE LINKS
echo
echo "** Creating symoblic links to libgptk source"
for i in `seq 0 $(( FILE_COUNT-1 ))`
do
  echo "Creating link ${LINK_NAME[$i]} -> $GPTK_SRC/${LINK_TO[$i]}"
  ln -s  $GPTK_SRC/${LINK_TO[$i]} ${LINK_NAME[$i]}
done

# config.h.in file
touch config.h.in
