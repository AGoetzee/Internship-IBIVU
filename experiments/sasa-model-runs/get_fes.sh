#!/bin/bash

echo
echo HILLS2FES processor v1.0
echo By Arthur Goetzee

ADDITIONAL_FLAGS=$1

if [ -z $ADDITIONAL_FLAGS ]
then
  echo "No additional flags given."
else
  echo "Additional flags given: $ADDITIONAL_FLAGS"
fi

for dir in */ 
do 

# Remove trailing /
dir=${dir%*/}

# Only works with numbers only
if [[ $dir =~ [0-9]+ ]]
then

  echo ****PROCESSING $dir******
  
  rm $i/fesdata -r
  mkdir $i/fesdata
  
  plumed sum_hills --hills $i/HILLS_$i --outfile $i/fesdata/fes.dat $ADDITIONAL_FLAGS
  plumed sum_hills --hills $i/HILLS_$i --stride 2000 --outfile $i/fesdata/fes_ $ADDITIONAL_FLAGS
  
  echo ****FINISHED $dir******
fi
done 

echo Finished processing hills. Exiting.
echo