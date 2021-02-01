"""
Some usefull function.
"""

import numpy as np
import synumses.one_dimension.parameters as parameters


def calc_p_density():
    """
    Calculates and returns an array of the hole density.
    """
    p_density = parameters.Nv*np.exp(parameters.q*( parameters.u[1::3] - parameters.u[0::3] - parameters.Chi - parameters.Eg)/(parameters.kB*parameters.T))
          
    return p_density



def calc_n_density():
    """
    Calculates and returns an array of the electron density.
    """
    n_density = parameters.Nc*np.exp(parameters.q*( parameters.u[0::3] - parameters.u[2::3] + parameters.Chi)/(parameters.kB*parameters.T))
    
    return n_density


def ohm_potential(C, Chi, Eg, Nc, Nv):
   """
   Calculates and returns the potential from the doping level, the electron affinity, the band gap, and the density of states.
   """

   from synumses.one_dimension.parameters import q, kB, Ut, T
   
   ni = np.sqrt(Nc*Nv)*np.exp(-((Eg)*q)/(2.*kB*T))
   Phi = -Chi - Eg/2. - 0.5*Ut*np.log(Nc/Nv) + Ut*np.arcsinh(C/(2.*ni))
    
   return(Phi)


def calc_ni_density():
   """
    Calculates and returns an array of the intrinsic carrier density.
    """
   ni_density = (
      np.sqrt(parameters.Nc*parameters.Nv)
      *
      np.exp(-((parameters.Eg)*parameters.q)
             /
             (2.*parameters.kB*parameters.T))
   )

   return ni_density

def calc_recombination():
   """
   Calculates  and returns an array of the recombination rate.
   
   """
   recombination = (
      parameters.Cau
      *
      (calc_p_density() * calc_n_density() -
       calc_ni_density()**2
      )
   )
    
   return recombination

def calc_generation(I_0, alpha):
   """
   Does nothing so far.
   """
   
   return None
    
