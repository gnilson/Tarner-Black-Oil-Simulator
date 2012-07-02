from pylab import *
from GAS import *
from OIL import *
from WATER import *
from ROCKPROPS import *
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
data["GORest"][0] = rsi = Rso(p, tres, pb, api, SGgas)


for i in range(1, len(pressures)):
    error = 1
    p = data["Pressure"][i]
    deltap = pi - p
    gor = data["GORest"][i-1] #new GOR guess
    rso = Rso(p, tres, pb, api, SGgas) #Solution GOR
    bo = Bo(p , tres, pb, rso, SGgas, api) * 5.615 #scf/stb
    bg = Bg(p, tres, SGgas) #rcuft/SCF
    et = bo - boi + bg * (rsi - rso) + boi * cer * deltap
    phiN = (bo - bg * rso) / et
    phiG = bg / et
    mug = Visg(p, tres, SGgas)
    muo = Visoil(p, tres, pb, SGgas, api)

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

        krg = Krg(sg, swi, sgc, sor, krg_sor)
        kro = Kro(so, swi, sor)
        if (krg < 0):
           krg = 0
        
        if (p<pi and p>pb):
            gorcalc = rsi
        elif (p<pb and sg<=sgc):
            gorcalc = rso
        else:
            gorcalc = rso + krg * muo \
                *(Bo(p, tres, pb, rso, SGgas, api)/5.615) / (kro \
                * mug * Bg(p, tres, SGgas))
        
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

figure(1, facecolor='w')

subplot(221)
plot(data["Pressure"], data["GORest"])
bottom, top = xlim()
xlim(top,bottom)
xlabel("Reservoir Pressure, psia")
ylabel("Producing GOR, SCF/STB")
grid(True)

subplot(222)
plot(data["Np"][1:], data["Pressure"][1:])
xlabel("Cumulative Oil Production, STB")
ylabel("Reservoir Pressure")
grid(True)

subplot(223)
plot(data["Pressure"][1:], data["PhiG"][1:])
bottom, top = xlim()
xlim(top,bottom)
xlabel("Reservoir Pressure, psia")
ylabel("Gas Saturation, fraction")
grid(True)

subplot(224)
plot(data["Pressure"][1:], data["Muo"][1:])
#bottom, top = xlim()
#xlim(top,bottom)
xlabel("Reservoir Pressure, psia")
ylabel("Oil Viscosity, Cp")
grid(True)

figure(2, facecolor='w')

subplot(221)
plot(data["Pressure"][1:], data["Mug"][1:])
xlabel("Reservoir Pressure, psia")
ylabel("Gas Viscosity, Cp")
grid(True)

subplot(222)
plot(data["Pressure"][1:], data["Muo"][1:])
xlabel("Reservoir Pressure, psia")
ylabel("Oil Viscosity, Cp")
grid(True)

subplot(223)
plot(data["Pressure"][1:], data["GORsol"][1:])
xlabel("Reservoir Pressure, psia")
ylabel("Solution Gas Oil Ratio, SCF/STB")
grid(True)


subplot(224)
plot(1-data["Sg"], data["krg"]/data["kro"])
xlabel("Liquid Saturation, Fraction")
ylabel("krg/kro")
grid(True)


show()
#print data["Muo"]
