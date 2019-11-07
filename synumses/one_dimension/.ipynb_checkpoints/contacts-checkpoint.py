import numpy as np

#from modells.parameters import *


def U_Ohm(C, Ec, Ev, Nc, Nv):

    from modells.parameters import q, kB, Ut, T
    #ni = np.sqrt(Nc*Nv)*np.exp(-((Ec-Ev)*parameters.q)/(2.*parameters.kB*parameters.T))
    #Phi = (Ec + Ev)/2. - 0.5*parameters.Ut*np.log(Nc/Nv) + parameters.Ut*np.arcsinh(C/(2.*ni))
    
    ni = np.sqrt(Nc*Nv*np.exp(-((Ec-Ev)*q)/(kB*T)))
    Phi = (Ec + Ev)/2. - 0.5*Ut*np.log(Nc/Nv) + Ut*np.arcsinh(C/(2.*ni))

    return(Phi)
