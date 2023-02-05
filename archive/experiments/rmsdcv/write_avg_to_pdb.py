import numpy as np
import matplotlib.pyplot as plt
from io import StringIO

# test = StringIO('0.000000 1.268718 1.022510 1.239088 0.847744 1.200293 0.894539 0.839608 0.878268 0.923797 1.058547 1.223668 1.208765 0.959895 0.886583 0.961046 0.800504 1.003612 0.935929 0.908011 0.829295 1.055866 1.024999 1.171880 0.964459 1.325895 1.042819 1.112009 0.816651 2.946817\n0.200000 1.217132 1.048282 1.244129 0.848938 1.189648 0.922272 0.860156 0.915951 0.955996 1.021675 1.236387 1.220793 0.998761 0.886463 1.039409 0.817188 0.921225 0.936838 0.955576 0.848323 1.049390 0.968498 1.185669 0.946163 1.314730 1.143355 1.141626 0.797883 2.994298')

data = np.loadtxt('rmsdcv/COLVAR_AVG').T
data = data[1:]
# data = np.loadtxt(test).T

avg = data.mean(axis=1)
sd = data.std(axis=1)
cs = ['c'+str(i+1) for i in range(29)]

plt.bar(cs,avg)
plt.errorbar(cs, avg, yerr=sd, fmt="o", color="r")
plt.xlabel('C-alpha atom')
plt.ylabel('Average distance (nm)')
# plt.show()
plt.savefig('rmsdcv\\avg_calpha_dist.png')

print('REMARK WEIGHT=1.0')
print(f'REMARK ARG={",".join(cs)}')
avg_list = []
avg = np.round(avg,decimals=2)
for i,(num_c,avg_c) in enumerate(zip(cs,avg)):
    avg_list.append(f'{num_c}={avg_c}')
print(f'REMARK {" ".join(avg_list)}')