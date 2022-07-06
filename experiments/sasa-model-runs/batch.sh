#!/bin/bash
#SBATCH --partition=gpu_shared
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=120:00:00

RESTART=false

while getopts 'rs:' OPTION; do

  case "$OPTION" in
    r)
			echo "Restarting"
			RESTART=true
			;;

    s)
			echo "The value of -s is $OPTARG"
			SCHEME=$OPTARG
			echo $SCHEME
      ;;

    ?)
      echo "Usage: $(basename $0) [-r] [-s argument]"
      exit 1
      ;;
  esac

done

CUDA_VISIBLE_DEVICES=0
CURRENT_DIR=`pwd`
TMPDIR=/scratch

date
echo ------Sourcing------
module load 2020
module load NCCL/2.7.8-gcccuda-2020a
module load intel/2020a

source /home/goetzee/amber20/amber.sh
source /home/goetzee/plumed/source.sh

echo All packages loaded


which plumed
which pmemd
module list

echo *******ENV VARS*********
printenv
echo *******ENV VARS*********

echo -------Copying-------
echo Copying to $TMPDIR

rsync -rv $SCHEME EM chain* *.pdb "$TMPDIR"
cd "$TMPDIR"

ls
date
pwd

echo -------Finished--------

if ! $RESTART
then 
  echo -------Heating---------
  echo Heating dir is $SCHEME/heat/
  srun pmemd.cuda -O -i $SCHEME/heat/heat.in -o $SCHEME/heat/heat_out.out -p chainsep_boxed.parm7 -c EM/EM_out.ncrst -r $SCHEME/heat/heat_out.ncrst -x $SCHEME/heat/heat_out.nc -inf $SCHEME/heat/heat_out.info
  echo Heating finished!
fi

echo ------Production-------
echo Production dir is $SCHEME/prod/
if $RESTART
then
  echo *****RESTARTING******
  srun pmemd.cuda -O -i $SCHEME/prod/restart.in -o $SCHEME/prod/prod_out_restart.out -p chainsep_boxed.parm7 -c $SCHEME/prod/restart.ncrst -r $SCHEME/prod/prod_out_restart.ncrst -x $SCHEME/prod/prod_out_restart.nc -inf $SCHEME/prod/prod_out_restart.info
else
  echo *****CLEAN START*****
  srun pmemd.cuda -O -i $SCHEME/prod/prod.in -o $SCHEME/prod/prod_out.out -p chainsep_boxed.parm7 -c $SCHEME/heat/heat_out.ncrst -r $SCHEME/prod/prod_out.ncrst -x $SCHEME/prod/prod_out.nc -inf $SCHEME/prod/prod_out.info
fi

echo Production finished.
echo -------Copying-------
ls
date

rsync -rv "$SCHEME" "$CURRENT_DIR"
echo -------Finished-------
date


