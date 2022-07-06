#!/bin/bash
# This script rsyncs a simulation.
echo
echo
echo "SyncSim 1.1"
echo "By Arthur Goetzee"

SCHEME=$1
FLAGS=$2
DEST=$3
JOBID=`squeue -u goetzee | awk -v scheme="$SCHEME" '$3 == scheme {print $1}'`
JOBNODE=`squeue -u goetzee | awk -v scheme="$SCHEME" '$3 == scheme {print $8}'`
RUNTIME=`squeue -u goetzee | awk -v scheme="$SCHEME" '$3 == scheme {print $6}'`
DEFAULTDEST="$HOME/simulations/sasa-model-runs/"
DEFAULTFLAGS="-vrun"


if ! [[ $JOBNODE =~ ^r[0-9]+n[0-9]+$ ]]
then
   echo "No JOBNODE found! Found $JOBNODE instead."
   echo "Exiting..." 
   exit 1
fi

if ! [[ $JOBID =~ ^[0-9]+$ ]]
then
   echo "No JOBID found!"
   echo "Exiting..."
   exit 1
fi


echo "SIMULATION NAME $SCHEME"
echo "Found JobID $JOBID"
echo "Found jobnode $JOBNODE"
echo "Running for $RUNTIME"
echo "Syncing from $JOBNODE:/scratch/slurm.$JOBID.0/scratch/$SCHEME" 
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
    FLAGS=$DEFAULTFLAGS
fi

echo "Flags given: $FLAGS"

# If no destination folder provided. Checks if $DEST exists.
if [ -z $DEST ]
then
	echo "No destination found. Using default: $DEFAULTDEST"
	echo
	echo "This will run the command:"
	echo "rsync $FLAGS $JOBNODE:/scratch/slurm.$JOBID.0/scratch/$SCHEME $DEFAULTDEST"
	run_check
	rsync $FLAGS $JOBNODE:/scratch/slurm.$JOBID.0/scratch/$SCHEME $DEFAULTDEST

else
	# If user gives destination folder
	echo "Destination folder found: $DEST"
	echo
	echo "This will run the command:"
	echo "rsync $FLAGS $JOBNODE:/scratch/slurm.$JOBID.0/scratch/$SCHEME $DEST"
	run_check
	rsync $FLAGS $JOBNODE:/scratch/slurm.$JOBID.0/scratch/$SCHEME $DEST
fi
