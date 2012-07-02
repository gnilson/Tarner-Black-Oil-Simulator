import GAS

import matplotlib.pyplot as plt
from numpy import *


ppr = arange(0, 8.1, 0.1)
tpr = arange(1.05, 3.05, 0.05)

Z = zeros((len(tpr), len(ppr)), dtype=float)

for i in range(len(tpr)):
    for j in range(len(ppr)):
        Z[i,j]= GAS.Zfactor(ppr[j], tpr[i])

fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111)
ax.grid(True)
for i in range(len(tpr)):
    ax.plot(ppr, Z[i],'k')

ax.set_xlabel("Pseudoreduced Pressure, Ppr")
ax.set_ylabel("Compressibility Factor, Z")
#figtext(.6, .3, "Pseudo Reduced Temperature")
plt.show()

