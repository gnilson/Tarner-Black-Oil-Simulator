from math import *

def Rso(p, t , pb , API , SGgas):
#Beggs and Standing Correlation
#UNITS=scf/STB
#T is reservoir temperature, F
#Rs is solution GOR in scf/STB
    Yg = 0.00091 * t - 0.0125 * API
    
    if (p > pb):
        return SGgas * (pb / (18 * 10 ** Yg)) ** 1.204
    else:
        return SGgas * (p / (18 * 10 ** Yg)) ** 1.204


def Visoil(p, t, pb, SGgas, API):

    #P -Pressure, psia
    #T -Temperature, F
    #pb-Bubble Point Pressure, psia
    #API- Oil API Gravity
    #Rsob-Solution gas oil ratio at bubble point
    #Rso-Solution Gas Oil Ratio
    #Function Returns as Oil viscosity
    
    myRso = Rso(p, t, pb, API, SGgas)
    myRsob = Rso(pb, t, pb, API, SGgas)

    A = 10.715 * ((myRsob + 100) ** (-0.515))
    B = 5.44 * ((myRsob + 150) ** (-0.338))
    Uob = A * (10 ** (10 ** (1.8653 - 0.025086 * API - 0.5644 * log10(t))) - 1) ** B
        

    if (p == pb):
        return Uob
    elif (p > pb):
        M = 2.6 * p ** 1.187 * (exp(-11.513 - ((8.98 * 10 ** (-5)) * p)))
        return Uob * ((p / float(pb)) ** M) 
    else:    
        C = 10.715 * (myRso + 100) ** (-0.515)
        D = 5.44 * (myRso + 150) ** (-0.338)
        Ucd = C * ((10 ** (10 ** (1.8653 - 0.025086 * API - 0.5644 * log10(t))) - 1) ** D)
        
    return Ucd


def Co(p, t, pb, sepGasGrav, apiOilGrav):
    rsob = Rso(pb, t, pb, apiOilGrav, sepGasGrav)
    

    if (p >= pb):
        A1 = -1433
        A2 = 5
        A3 = 17.2
        A4 = -1180
        A5 = 12.61
        A6 = 10 ** 5
        
        return (A1 + A2 * rsob + A3 * t + A4 * sepGasGrav + A5 * apiOilGrav) / (A6 * p)
    
    else:
        #This is the Villena and Lanzi correlation from Craft&Hawkins pg 39
        return exp(-0.664 - 1.43 * log(p) - 0.395 * log(pb) + 0.39 * log(t+459.67) + 0.455 * log(rsob) + 0.262 * log(apiOilGrav))

    
def Bo(p, t, pb, rso, sepGasGrav, apiOilGrav): 
    sgo = (141.5 / (131.5 + apiOilGrav))
    f = rso * (sepGasGrav / sgo) ** 0.5 + 1.25 * t
    bob = 0.972 + 0.000147 * f ** 1.175
      
    if (p > pb):
        return bob * exp(Co(p, t, pb, sepGasGrav, apiOilGrav) * (pb - p))
    else:
        return bob


def Pb_almarhoun(rsp, SGgas, sgSTO, t):

    """Al-Marhoun "PVT Correlations for Middle East Crude Oils" from McCain pg 519
    Note that the equation gives better results when the separator gas-oil ratio is used rather than
    the total gas-oil ratio."""
    
    temp = t + 459.67
    
    a1 = 5.38088 * 10 ^ -3
    a2 = 0.715082
    a3 = -1.87784
    a4 = 3.1437
    a5 = 1.32657
    
    return a1 * (rsp ** a2) * (SGgas ** a3) * (sgSTO ** a4) * (temp ** a5)

def Pb_standing(rs, SGgas, APIoil, t):

    """Standing bubble point correlation
    Rs should include separator gas and stock-tank gas
    T in deg F"""
    
    return 18.2 * ((rs / SGgas) ** 0.83 * 10 ** (0.00091 * t - 0.0125 * APIoil) - 1.4)
    
