import numpy as np
import matplotlib.pyplot as plt

data_temp = np.loadtxt('summary.TEMP')

plt.plot(data_temp[:,0],data_temp[:,1])
plt.xlabel('Time (ps)')
plt.ylabel('Temperature (K)')
plt.title('Temperature over time')
plt.savefig('003_temp.png')
plt.clf()

data_energy1 = np.loadtxt('summary.ETOT')
data_energy2 = np.loadtxt('summary.EKTOT')
data_energy3 = np.loadtxt('summary.EPTOT')

plt.plot(data_energy1[:,0],data_energy1[:,1], label='Total energy')
plt.plot(data_energy2[:,0],data_energy2[:,1], label='Kinetic energy')
plt.plot(data_energy3[:,0],data_energy3[:,1], label='Potential energy')
plt.legend()
plt.xlabel('Time (ps)')
plt.ylabel('Energy (kcal/mol)')
plt.title('energy over time')
plt.savefig('003_energy.png')
plt.clf()
