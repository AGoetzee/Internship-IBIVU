#!/bin/bash

# Generate topology
printf '6\n1' | gmx pdb2gmx -ignh -f 6nzn.pdb

# Put the protein in a box
gmx editconf -f conf.gro -bt dodecahedron -d 1 -o box.gro

# Put water in the box
gmx solvate -cs -cp box.gro -p topol.top -o solve.gro

# Energy minimization
gmx grompp -f em.mdp -c solve.gro -p topol.top -o emw.tpr -maxwarn 1

# Run EM
gmx mdrun -deffnm emw -v

# Add ions
gmx grompp -f em.mdp -c emw.gro -p topol.top -o ion.tpr -maxwarn 1

printf '13' | gmx genion -s ion.tpr -p topol.top -neutral -conc 0.1 -o ion.gro

# Energy minimize again
gmx grompp -f em.mdp -c ion.gro -p topol.top -o em.tpr

gmx mdrun -deffnm em -v

# Run position restraints
gmx grompp -f posre.mdp -c em.gro -p topol.top -o posre.tpr -r em.gro

gmx mdrun -deffnm posre -v

# Actual MD run
gmx grompp -f md.mdp -c posre.gro -p topol.top -o md.tpr

gmx mdrun -deffnm md -v