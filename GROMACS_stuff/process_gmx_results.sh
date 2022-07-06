#!/bin/bash

# GROMACS Result processor
# Arthur Goetzee 9/10/21
# Usage: Exectute in experiment folder OR execute in any folder and supply experiment path as first arg.

module load 2020
module load GROMACS/2020.3-intelcuda-2020a

# Check if folder is given. If it is, go that folder.
if [ -n "$1" ]; then
    echo "Folder given. Going to exp folder..."
	cd $1
else
	echo "No folder given. Executing in the current folder."
fi

echo Processing results...

# Get current exp name
SIM=$(basename `pwd`)
echo Experiment $SIM

# Set output dir
RESULTS="../../results/$SIM"

# Create output dir
# Note: May throw an error if dir already exists. Can be ignored
mkdir $RESULTS

# Get all energy stats
printf '1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n21\n22\n23\n24\n25\n26\n27\n28\n29\n30\n31\n32\n33\n34\n35\n36\n37\n38\n39\n40\n41\n42\n43\n44\n45\n46\n47\n48\n49\n50\n51\n\n' | gmx energy -f $SIM.edr -o $RESULTS/energy.xvg

# Get internal hydrogen bonds
printf '1\n1' | gmx hbond -f $SIM.xtc -s $SIM.tpr -num $RESULTS/hbonds-int.xvg

# Get solvent-protein hbonds
printf '1\n13' | gmx hbond -f $SIM.xtc -s $SIM.tpr -num $RESULTS/hbonds-ext.xvg

# Get RMSD to native structure
printf '3' | gmx rms -f $SIM.xtc -s $SIM.tpr -o $RESULTS/rmsd.xvg

# Get SASA
printf '1\n' | gmx sasa -f $SIM.xtc -s $SIM.tpr -o $RESULTS/sasa.xvg

# Copy relevant files to results folder
cp $SIM-batch.sh $RESULTS

echo Results of $SIM processed. Saved to $RESULTS