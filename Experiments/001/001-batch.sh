#!/bin/bash
#SBATCH -N 2 --tasks-per-node=16
#SBATCH -t 6:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=arthurgoetzee@outlook.com

module load 2020
module load GROMACS/2020.3-intelcuda-2020a

# Process flags
# -p: path to structure file
# -s: Experiment number
while getopts p:s: flag
do
    case "${flag}" in
        p) PDBFILE=${OPTARG};;
        s) SIM=${OPTARG};;
    esac
done
echo "Starting Experiment $SIM"
echo "Using structure file $PDBFILE"


### Generate topology ###
# Atomistic forcefield and explicit water
printf '6\n1' | gmx pdb2gmx -ignh -f $PDBFILE

### Put the protein in a box ###
# Distance 1
# Box is a dodecahedron
gmx editconf -f conf.gro -bt dodecahedron -d 1 -o box.gro

### Put water in the box ##
gmx solvate -cs -cp box.gro -p topol.top -o solve.gro

### Add ions ###
# Check later if this is necessary
gmx grompp -f ions.mdp -c solve.gro -p topol.top -o ion.tpr -maxwarn 1

printf '13' | gmx genion -s ion.tpr -p topol.top -neutral -conc 0.1 -o ion.gro

### Energy minimize ###
gmx grompp -f em.mdp -c ion.gro -p topol.top -o em.tpr

gmx mdrun -deffnm em -v

### Equilibration ###
# NVT
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

gmx mdrun -deffnm nvt -v

# NPT
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

gmx mdrun -deffnm npt

# Actual MD run
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o $SIM.tpr

gmx mdrun -deffnm $SIM

# Process results
../../scripts/process_gmx_results.sh