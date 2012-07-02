from math import *

A1 = 0.3265
A2 = -1.07
A3 = -0.5339
A4 = 0.01569
A5 = -0.05165
A6 = 0.5475
A7 = -0.7361
A8 = 0.1844
A9 = 0.1056
A10 = 0.6134
A11 = 0.721


def Bg(p, t, SGgas):
    #cu. ft./SCF
    Tpr = (t + 459.67)/Tpc(SGgas)
    Ppr = p / Ppc(SGgas)
    return 0.02829 * (Zfactor(Ppr, Tpr) * (t + 459.67))/p

def Tpc(SGgas):
    return 169.2 + (349.5 * SGgas) - (74 * SGgas**2)

def Tpc_condensate(SGgas):
    """Sutton 2007 SPE Reservoir Evaluation & Engineering 2008 vol 10 issue 3"""
    return 164.3 + (357.7 * SGgas) - (67.7 * (SGgas ** 2))

def Ppc(SGgas):
    return 756.8 - (131 * SGgas) - (3.6 * SGgas**2)

def Ppc_condensate(SGgas):
    """Sutton 2007 SPE Reservoir Evaluation & Engineering 2008 vol 10 issue 3"""
    return 744 - (125.4 * SGgas) + (5.9 * (SGgas ** 2))

def Ppr(p, SGgas):
    return p/Ppc(SGgas)

def Tpr(t, SGgas):
    return (t + 459.67)/Tpc(SGGas)

def Cg(p, t, SGgas):
    Tpr = (t + 459.67) / Tpc(SGgas)
    Ppr = p / Ppc(SGgas)
    z = Zfactor(Ppr, Tpr)
    rhor = 0.27 * Ppr / (z * Tpr)
    
    dZdRHO = A1 + (A2 / Tpr) + (A3 / Tpr ** 3) + (A4 / Tpr ** 4) \
        + (A5 / Tpr ** 5) + 2 * rhor * (A6 + (A7 / Tpr) + (A8 / Tpr ** 2)) \
        - 5 * rhor ** 4 * A9 * ((A7 / Tpr) + (A8 / Tpr ** 2)) \
        + ((2 * A10 * rhor) / Tpr ** 3) * (1 + A11 * rhor ** 2 \
        - A11 ** 2 * rhor ** 4) * exp(-A11 * rhor ** 2)
             
    cr = (1 / Ppr) - (0.27 / (z ** 2 * Tpr)) * (dZdRHO / (1 + (rhor / z) * dZdRHO))
    return cr / Ppc(SGgas)

def Zfactor(ppr, tpr):
    ztolerance = 0.0001
        
    c1 = A1 + A2 / tpr + A3 / tpr ** 3 + A4 / tpr ** 4 + A5 / tpr ** 5
    c2 = A6 + A7 / tpr + A8 / tpr ** 2
    c3 = A9 * (A7 / tpr + A8 / tpr ** 2)
    
    z = 0.9 #Z guess
    rhor = 0.27 * ppr / (z * tpr)
    c4 = A10 * (1 + A11 * rhor ** 2) \
        * (rhor ** 2 / tpr ** 3) \
        * exp(-A11 * rhor ** 2)
    f = z - (1 + c1 * rhor + c2 * rhor ** 2 - c3 * rhor ** 5 + c4)
        
    while (abs(f) > ztolerance):
        df = 1 + c1 * rhor / z + 2 * c2 * rhor ** 2 / z - 5 * c3 * rhor ** 5 / z \
            + ((2 * A10 * rhor ** 2) / (tpr ** 3 * z)) \
            *(1 + A11 * rhor ** 2 - (A11 * rhor ** 2) ** 2) \
            * exp(-A11 * rhor ** 2)
        z = z - f / df
        rhor = 0.27 * ppr / (z * tpr)
        c4 = A10 * (1 + A11 * rhor ** 2) * \
             (rhor ** 2 / tpr ** 3) * \
             exp(-A11 * rhor ** 2)
        f = z - (1 + c1 * rhor + c2 * rhor ** 2 - c3 * rhor ** 5 + c4)
             
    return z

def Visg(p, t, SGgas):
    myT = t + 459.67
    M = 29 * SGgas
    ppr = p / Ppc(SGgas)
    tpr = (t + 459.67) / Tpc(SGgas)
    rhog = 1.4935 * (10 ** -3) * p * M / (Zfactor(ppr, tpr) * myT)
    
    A = ((9.379 + 0.01607 * M) * myT ** 1.5) / (209.2 + 19.26 * M + myT)
    B = 3.448 + (986.4 / myT) + (0.01009 * M)
    C = 2.447 - (0.2224 * B)
    
    return A * (10 ** -4) * exp(B * (rhog ** C))

def GasDensity(p, t, SGgas):
    ppr = p / Ppc(SGgas)
    tpr = (t + 459.67) / Tpc(SGgas)
    z = Zfactor(ppr, tpr)
    Ma = (SGgas * 28.97)
    B = (t + 459.67)
    C = (453.59237 / (30.48 ^ 3))
 
    return ((p * Ma) / (10.732 * B * z)) * C
    

def StaticBottomHolePressure(Ptop, Ttop, Tbott, SGgas, depth):
    """Calculation of pressure change due to column of gas of a given height"""

    error = 20
    pguess = Ptop + 0.07 * depth
    ztop = Zfactor(Ppr(Ptop, SGgas), Tpr(Ttop, SGgas))
    tave = (Ttop + Tbott) / 2 + 459.67

    while (error > 0.01):
        zguess = Zfactor(Ppr(pguess, SGgas), Tpr(Tbott, SGgas))
        zave = (zguess + ztop) / 2
        pcalc = Ptop * Exp((0.01875 * SGgas * depth) / (zave * (tave + 459.67)))
        error = Abs(pcalc - pguess)
        pguess = pcalc

    return pcalc


