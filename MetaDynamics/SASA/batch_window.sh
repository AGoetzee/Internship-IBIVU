#!/bin/bash
#SBATCH --partition=binf,defq
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00

# This is the BAZIS version

RESTART=false

while getopts 'ret:c:s:' OPTION; do

  case "$OPTION" in
    r)
			echo "Restarting"
			RESTART=true
			;;

    e)
            echo "Energy minimizing"
            EM=true
            ;;
            
    t)
            echo "The topology is $OPTARG"
            TOPOL=$OPTARG
            ;;
    c)
            echo "The starting config is $OPTARG"
            CONFIG=$OPTARG
            ;;
            
    s)
			SCHEME=$OPTARG
			echo "The scheme is $SCHEME"
      ;;

    :)
      echo "$0: Must supply an argument to -$OPTARG." >&2
      exit 1
      ;;

    ?)
      echo "Usage: $(basename $0) [-r] [-e] [-s simulation scheme] [-t topology (.parm7)] [-c coordinate file (.rst7)]"
      exit 1
      ;;
  esac

done

CUDA_VISIBLE_DEVICES=0
CURRENT_DIR=`pwd`

date
echo ------Sourcing------

module load amber/20
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

rsync -rv "$SCHEME" src "$TMPDIR"
cd "$TMPDIR"

ls
date
pwd

echo -------Finished--------

if $EM
then
  echo ----------EM------------
  echo EM dir is $SCHEME/EM
  srun pmemd.cuda -O -i $SCHEME/EM/EM.in -o $SCHEME/EM/EM_out.out -p $TOPOL -c $CONFIG -r $SCHEME/EM/EM_out.ncrst -inf $SCHEME/EM/EM_out.info
  echo EM finished!
fi
  
if ! $RESTART
then 
  echo 
  echo Heating dir is $SCHEME/heat/
  srun pmemd.cuda -O -i $SCHEME/heat/heat.in -o $SCHEME/heat/heat_out.out -p $TOPOL -c $SCHEME/EM/EM_out.ncrst -r $SCHEME/heat/heat_out.ncrst -x $SCHEME/heat/heat_out.nc -inf $SCHEME/heat/heat_out.info
  echo Heating finished!
fi

echo ------Production-------
echo Production dir is $SCHEME/prod/
if $RESTART
then
  echo *****RESTARTING******
  srun pmemd.cuda -O -i $SCHEME/prod/restart.in -o $SCHEME/prod/prod_out_restart.out -p $TOPOL -c $SCHEME/prod/restart.ncrst -r $SCHEME/prod/prod_out_restart.ncrst -x $SCHEME/prod/prod_out_restart.nc -inf $SCHEME/prod/prod_out_restart.info
else
  echo *****CLEAN START*****
  srun pmemd.cuda -O -i $SCHEME/prod/prod.in -o $SCHEME/prod/prod_out.out -p $TOPOL -c $SCHEME/heat/heat_out.ncrst -r $SCHEME/prod/prod_out.ncrst -x $SCHEME/prod/prod_out.nc -inf $SCHEME/prod/prod_out.info
fi

echo Production finished.
echo -------Copying-------
ls
date

rsync -rv "$SCHEME" "$CURRENT_DIR"
echo -------Finished-------
date


