__author__= 'Arthur Goetzee'

from pytraj import load, write_traj # For Amber trajectories
import numpy as np # Data analysis
from os import walk, getcwd, path
from subprocess import call

def check_wd(wd):
    print(f'Current dir: {wd}')
    if not ((wd == '/scistor/informatica/age206/simlab/final-script') or (wd == '/home/goetzee/simulations/sasa-model-runs')):
        print('You are in the wrong working directory!\nExecuting this in the wrong directory can lead to unpredictable behaviour.\n')
        print('Please navigate to the following:')
        print('BAZIS:/scistor/informatica/age206/simlab/final-script')
        print('LISA: /home/goetzee/simulations/sasa-model-runs\n')
        print('Exiting.')
        exit()

print("\n\tgenerate_restart.py\nGenerate restart frames including velocities.\nMade by Arthur Goetzee, ver 1.1\n")

check_wd(getcwd())

sim_dirs = []
print('Going through dirs....')
for root,dirs,files in walk('.'):
    for dir in dirs:
        if dir.isnumeric():
           print(f'Found {dir}')
           sim_dirs.append(dir)


for sim in sim_dirs:
    print(f'Processing {sim}...')
    
    if not path.exists(f'{sim}/prod/prod_out.nc'):
        print(f'Trajectory does not exist, Skipping...')
        continue
    
    traj = load(f'{sim}/prod/prod_out.nc',"chainsep_boxed.parm7")
    write_traj(f"{sim}/prod/restart.ncrst",traj,frame_indices=[traj.n_frames-1],time=True,overwrite=True,velocity=True)
    call(["mv",f"{sim}/prod/restart.ncrst.1",f"{sim}/prod/restart.ncrst"])
    print(f'File written to {sim}/prod/restart.ncrst')

print('Finished!\n')
