import pytraj as pt
import numpy as np
import matplotlib.pyplot as plt

traj = pt.load('002_Pin.nc', top='002.parm7')
rmsd_data = pt.rmsd(traj, ref='002_Min.ncrst')

plt.plot(rmsd_data)
plt.xlabel('Frame')
plt.ylabel('RMSD')
plt.title('RMSD per frame')
plt.savefig('002_rmsd.png')
