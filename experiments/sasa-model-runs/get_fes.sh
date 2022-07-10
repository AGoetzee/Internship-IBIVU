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
  
  rm $dir/fesdata -r
  mkdir $dir/fesdata
  
  plumed sum_hills --hills $dir/HILLS_$dir --outfile $dir/fesdata/fes.dat $ADDITIONAL_FLAGS
  plumed sum_hills --hills $dir/HILLS_$dir --stride 2000 --outfile $dir/fesdata/fes_ $ADDITIONAL_FLAGS
  
  echo ****FINISHED $dir******
fi
done 

echo Finished processing hills. Exiting.
echo