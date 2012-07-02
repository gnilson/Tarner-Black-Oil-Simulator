import matplotlib.pyplot as plt
import GAS
import OIL
import WATER
import ROCKPROPS


from numpy import *
from PARAMETERS import *

pressures = []
for p in range(pi, pab-pstep, -pstep):
    pressures.append(p)
    if (p>pb and (p-pstep)<pb):
        pressures.append(pb)
    
data = {"Pressure":array(pressures),
        "GORsol":zeros(len(pressures), 'f'),
        "GORest":zeros(len(pressures), 'f'),
        "GORave":zeros(len(pressures), 'f'),
        "PhiN":zeros(len(pressures), 'f'),
        "PhiG":zeros(len(pressures), 'f'),
        "deltaNp":zeros(len(pressures), 'f'),
        "Np":zeros(len(pressures), 'f'),
        "ORF":zeros(len(pressures), 'f'),
        "So":zeros(len(pressures), 'f'),
        "Sg":zeros(len(pressures), 'f'),
        "kro":zeros(len(pressures), 'f'),
        "krg":zeros(len(pressures), 'f'),
        "deltaGp":zeros(len(pressures), 'f'),
        "Gp":zeros(len(pressures), 'f'),
        "Muo":zeros(len(pressures), 'f'),
        "Mug":zeros(len(pressures), 'f')}


p = pi
cer = (cw * swi + cf)/(1-swi)
boi = boi * 5.615 #scf/stb
data["GORest"][0] = rsi = OIL.Rso(p, tres, pb, api, SGgas)


for i in range(1, len(pressures)):
    error = 1
    p = data["Pressure"][i]
    deltap = pi - p
    gor = data["GORest"][i-1] #new GOR guess
    rso = OIL.Rso(p, tres, pb, api, SGgas) #Solution GOR
    bo = OIL.Bo(p , tres, pb, rso, SGgas, api) * 5.615 #scf/stb
    bg = GAS.Bg(p, tres, SGgas) #rcuft/SCF
    et = bo - boi + bg * (rsi - rso) + boi * cer * deltap
    phiN = (bo - bg * rso) / et
    phiG = bg / et
    mug = GAS.Visg(p, tres, SGgas)
    muo = OIL.Visoil(p, tres, pb, SGgas, api)

    while(error>tolerance):
        gorave = (gor + data["GORest"][i-1])/2
        deltaNp = (N - data["Np"][i-1]*phiN - data["Gp"][i-1]*phiG) \
            / (phiN + gorave*phiG)
        Np = data["Np"][i-1] + deltaNp
        deltaGp = gorave * deltaNp
        Gp = data["Gp"][i-1] + deltaGp
        ORF = Np / N
        so = (1 - Np/N)*(1-swi)*bo/(boi*(1-cf*deltap))
        sw = swi * (1 + cw*deltap)/(1-cf*deltap)
        
        #if (p>=pb):
        #    sg = 0
        #else:
        sg = 1 - so - sw

        krg = ROCKPROPS.Krg(sg, swi, sgc, sor, krg_sor)
        kro = ROCKPROPS.Kro(so, swi, sor)

        if (krg < 0):
           krg = 0
        if (p<pi and p>pb):
            gorcalc = rsi
        elif (p<pb and sg<=sgc):
            gorcalc = rso
        else:
            gorcalc = rso + krg * muo \
                *(OIL.Bo(p, tres, pb, rso, SGgas, api)/5.615) / (kro \
                * mug * GAS.Bg(p, tres, SGgas))
        
        error = abs(gorcalc - gor)/gorcalc
        gor = gorcalc

    data["GORsol"][i] = rso
    data["GORest"][i] = gor
    data["GORave"][i] = (data["GORest"][i] + data["GORest"][i-1])/2
    data["PhiN"][i] = phiN
    data["PhiG"][i] = phiG
    data["deltaNp"][i] = deltaNp
    data["Np"][i] = Np
    data["ORF"][i] = ORF
    data["So"][i] = so
    data["Sg"][i] = sg
    data["kro"][i] = kro
    data["krg"][i] = krg
    data["Gp"][i] = Gp
    data["Muo"][i] = muo
    data["Mug"][i] = mug

fig1 = plt.figure(1, facecolor='w')

ax = fig1.add_subplot(221)
ax.plot(data["Pressure"], data["GORest"])
#bottom, top = plt.xlim()
#plt.xlim(top,bottom)
ax.set_xlabel("Reservoir Pressure, psia")
ax.set_ylabel("Producing GOR, SCF/STB")
ax.grid(True)

ax = fig1.add_subplot(222)
ax.plot(data["Np"][1:], data["Pressure"][1:])
ax.set_xlabel("Cumulative Oil Production, STB")
ax.set_ylabel("Reservoir Pressure")
ax.grid(True)

ax = fig1.add_subplot(223)
ax.plot(data["Pressure"][1:], data["PhiG"][1:])
#bottom, top = xlim()
#ax.xlim(top,bottom)
ax.set_xlabel("Reservoir Pressure, psia")
ax.set_ylabel("Gas Saturation, fraction")
ax.grid(True)

ax = fig1.add_subplot(224)
ax.plot(data["Pressure"][1:], data["Muo"][1:])
#bottom, top = xlim()
#xlim(top,bottom)
ax.set_xlabel("Reservoir Pressure, psia")
ax.set_ylabel("Oil Viscosity, Cp")
ax.grid(True)

fig2 = plt.figure(2, facecolor='w')

ax = fig2.add_subplot(221)
ax.plot(data["Pressure"][1:], data["Mug"][1:])
ax.set_xlabel("Reservoir Pressure, psia")
ax.set_ylabel("Gas Viscosity, Cp")
ax.grid(True)

ax = fig2.add_subplot(222)
ax.plot(data["Pressure"][1:], data["Muo"][1:])
ax.set_xlabel("Reservoir Pressure, psia")
ax.set_ylabel("Oil Viscosity, Cp")
ax.grid(True)

ax = fig2.add_subplot(223)
ax.plot(data["Pressure"][1:], data["GORsol"][1:])
ax.set_xlabel("Reservoir Pressure, psia")
ax.set_ylabel("Solution Gas Oil Ratio, SCF/STB")
ax.grid(True)

ax = fig2.add_subplot(224)
ax.plot(1-data["Sg"], data["krg"]/data["kro"])
ax.set_xlabel("Liquid Saturation, Fraction")
ax.set_ylabel("krg/kro")
ax.grid(True)


plt.show()
#print data["Muo"]
