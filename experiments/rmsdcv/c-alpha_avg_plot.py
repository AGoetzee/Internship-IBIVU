import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('rmsdcv\\run2\\COLVAR_run2').T
data = data[1:28]
# data = np.loadtxt(test).T

avg = data.mean(axis=1)
sd = data.std(axis=1)
cs = ['c'+str(i+1) for i in range(27)]

plt.bar(cs,avg)
plt.errorbar(cs, avg, yerr=sd, fmt="o", color="r")
plt.xlabel('C-alpha atom')
plt.ylabel('Average distance (nm)')
# plt.show()
plt.savefig('rmsdcv\\avg_calpha_dist_run2.png')

print('REMARK WEIGHT=1.0')
print(f'REMARK ARG={",".join(cs)}')
avg_list = []
avg = np.round(avg,decimals=2)
for i,(num_c,avg_c) in enumerate(zip(cs,avg)):
    avg_list.append(f'{num_c}={avg_c}')
print(f'REMARK {" ".join(avg_list)}')