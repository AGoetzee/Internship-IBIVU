import pytraj as pt
import numpy as np
import matplotlib.pyplot as plt

traj = pt.load('003_Pin.nc', top='003.parm7')

rmsd_data = pt.rmsd(traj, ref=0)

plt.plot(rmsd_data)
plt.xlabel('Frame')
plt.ylabel('RMSD')
plt.title('RMSD per frame')
plt.savefig('003_rmsd.png')
