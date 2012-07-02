from math import *

def Kro(So , Swi , Sorg): 
    return 0.98372 * (So / (1 - Swi)) ** 4 * ((So - Sorg) / (1 - Swi - Sorg)) ** 2
    

def Krg(Sg , Swi , Sgc , Sorg , KrgAtSorg): 
    if (Sg==0):
        return 0
    else:
        return 1.1072 * ((Sg - Sgc) / (1 - Swi)) ** 2 * KrgAtSorg + 2.7794 \
            * ((Sorg * (Sg - Sgc)) / (1 - Swi)) * KrgAtSorg

def Cf(phi):
    #Formation Compressibility Correlation
    
    return (97.32 * 10 ** -6) / (1 + 55.8721 * phi) ** 1.42859
        

