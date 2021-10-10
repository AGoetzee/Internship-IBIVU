Each experiment has 3 digits, given in the order that they appear.

The experiment folder can be found in ~/Experiments/<Exp digit code>

Each experiment folder contains:
- A word file with a short description
- .mdp files
- .top files
- .gro files
- .ndx files
- Bash script used (if it was modified; check labjournal)

The results can be found in ~/Results/<Exp digit code>

Each result folder contains the following, by default: (Generated from automated bash scripts. see ~/Scripts)
- MobaXTerm terminal output (only contains the execution up until the batch script execution on LISA)
- energy.xvg containing -all- energy statistics vs time. Obtained through gmx energy.
- sasa.xvg containing SASA of the protein vs time. Obtained through gmx sasa.
- rms.xvg containg RMSD to the reference structure vs time.  Obtained through gmx rms.
- hbond-int.xvg containing total internal protein hbonds vs time. Obtained through gmx hbond.
- hbond-sol.xvg containing total external protein hbonds with the solvent vs time. Obtained through gmx hbond.

Other more specific results should be described in the labjournal and in the <exp digit code> folder.
