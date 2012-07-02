from GAS import *
from pylab import *

ppr = arange(0, 8.1, 0.1)
tpr = arange(1.05, 3.05, 0.05)

Z = zeros((len(tpr), len(ppr)), dtype=float)

for i in range(len(tpr)):
    for j in range(len(ppr)):
        Z[i,j]= Zfactor(ppr[j], tpr[i])

figure(1, facecolor='w')
grid(True)
for i in range(len(tpr)):
    plot(ppr, Z[i],'k')

xlabel("Pseudoreduced Pressure, Ppr")
ylabel("Compressibility Factor, Z")
#figtext(.6, .3, "Pseudo Reduced Temperature")
show()
print Z[3]

