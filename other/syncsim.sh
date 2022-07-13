#!/bin/bash
# This script rsyncs a simulation.
echo
echo
echo "SyncSim 1.1"
echo "By Arthur Goetzee"

SCHEME8=${1::8}
SCHEME=$1
FLAGS=$2
DEST=$3
JOBID=`squeue -u $USER | awk -v scheme="$SCHEME8" '$3 == scheme {print $1}'`
JOBNODE=`squeue -u $USER | awk -v scheme="$SCHEME8" '$3 == scheme {print $8}'`
RUNTIME=`squeue -u $USER | awk -v scheme="$SCHEME8" '$3 == scheme {print $6}'`
DEFAULTFLAGS="-vrun"


# Check which HPC the script runs on
if [[ `hostname` =~ surfsara ]]
then
    echo "LISA Node detected."
    DEFAULTDEST="$HOME/simulations/sasa-model-runs/"
    SOURCE="$JOBNODE:/scratch/slurm.$JOBID.0/scratch/$SCHEME"

else
    echo "BAZIS Node detected."
    DEFAULTDEST="$HOME/simlab/final-script/"
    SOURCE="$JOBNODE:/scratch/$JOBID/$SCHEME"
fi

if ! [[ $JOBNODE =~ ^r[0-9]+n[0-9]+$ ]] && ! [[ $JOBNODE =~ node[0-9]+ ]]
then
   echo "No JOBNODE found! Found '$JOBNODE' instead."
   echo "Exiting..."
   exit 1
fi

if ! [[ $JOBID =~ ^[0-9]+$ ]]
then
   echo "No JOBID found! Found '$JOBID' instead."
   echo "Exiting..."
   exit 1
fi


echo "SIMULATION NAME $SCHEME"
echo "Found JobID $JOBID"
echo "Found jobnode $JOBNODE"
echo "Running for $RUNTIME"
echo "Syncing from $SOURCE"
echo

run_check () {
read -p "Are you sure? " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    exit 1
fi
}

if [ -z $FLAGS ]
then
    echo "No flags given, uing default: $DEFAULTFLAGS"
    echo "WARNING: This instance will NOT sync!"
    FLAGS=$DEFAULTFLAGS
else
    echo "Flags given: $FLAGS"
fi


# If no destination folder provided. Checks if $DEST exists.
if [ -z $DEST ]
then
	echo "No destination found. Using default: $DEFAULTDEST"
	echo
	echo "This will run the command:"
	echo "rsync $FLAGS $SOURCE $DEFAULTDEST"
	run_check
	rsync $FLAGS $SOURCE $DEFAULTDEST

else
	# If user gives destination folder
	echo "Destination folder found: $DEST"
	echo
	echo "This will run the command:"
	echo "rsync $FLAGS $SOURCE $DEST"
	run_check
	rsync $FLAGS $SOURCE $DEST
fi
