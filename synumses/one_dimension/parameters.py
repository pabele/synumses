"""
Within this module all global variables and arrays are defined.
It also defines the default device, a silicon pn junction.

.. csv-table:: **global variables**
   :header: "global variable", "meaning", "default value"
   :widths: 10, 25, 10
   
   "T", "temperature", "300"
   "q", "elementary charge :math:`q_\mathrm{e}`", "1.6E-19"
   "kB", "Boltzman constant :math:`k_\mathrm{B}`",  "1.38E-23"
   "Ut", "temperature voltage :math:`U_\mathrm{T}`", "kB*T/q"
   "Epsilon_0", "vacuum permittivity :math:`\epsilon_0`", "8.85E-12"
   "Epsilon_r", "relative permittivity of silicon :math:`\epsilon_\mathrm{r}`",  "11.6"
   "n", "number of grid points", "400"
   "lx", "length of device", "400E-9"
   "bernoulli_limit", "for values less then **bernoulli_limit** the bernoulli function is approximated by a polynom", "1E1"

.. csv-table:: **global arrays**
   :header: "global array", "meaning", "default value"
   :widths: 10, 25, 10
   
   "pos_x", "position of grid point", "np.linspace(0,lx,n)"
  
   "Chi", "electron affinity", "np.full(n, 4.05)"
   "Eg", "band gap", "np.full(n, 1.12)"


   "Nc", "effective density of states in conduction band", "np.full(n, 2.81E25)"
   "Nv", "effective density of states in valence band", "np.full(n, 1.83E25)"


   "Epsilon", "permittivity of the material", "np.full(n, Epsilon_r * Epsilon_0)"
   "mu_p", "hole mobility :math:`\mu_\mathrm{p}`", "np.full(n, 0.045)"
   "mu_n", "electron mobility :math:`\mu_\mathrm{n}`", "np.full(n, 0.14)"
   "C", "doping", "np.zeros(n), -1E24/1E24"
   "Cau", "coefficient for recombination: q*(Cau*(n*p-ni2)-generation)*dx ", "np.full(n, 0)"
   "generation", "generation: q*(Cau*(n*p-ni2)-generation)*dx", "np.full(n, 0.0)"

   "b", ":math:`f(x)`", "np.zeros(3*n)"
   "A", "jacobian matrix",  "sparse.lil_matrix((3*n, 3*n))"
   "x", ":math:`\delta x`",  "np.zeros(3*n)"

"""


import numpy as np

import scipy.sparse as sparse 


# Global Variables

# #####################
#       Constants
# #####################
T = 300.
q = 1.6E-19

kB = 1.38E-23
Ut = kB*T/q

Epsilon_0 = 8.85E-12
Epsilon_r = 11.6

# ############
#    Mesh
# ############
n  = 400
lx = 400E-9


# *********************************************
#  Limit for Bernoulli function epproximation
# *********************************************

bernoulli_limit = 1E-4 # Default 1E-4


def init_geometry():
        """
        This functions sets grid size (**dx**) and the position of the gris points (**pos_x**).
        This function must be executed after changing the length (**lx**) or the number of grid points (**n**).    
        """
        

        global dx
        dx = lx/n

        global pos_x
        pos_x = np.linspace(0,lx,n)

def init_parameters():
        """
        This function defines the material parameters for silicon.
        After executing this function the parameters can be altered.
        This function must be executed after **init_geometry()**
        """

        global dt
        dt = 1E-12
        
        global Chi
        Chi = np.full(n, 4.05)

        global Eg
        Eg = np.full(n, 1.12)
        
        global Nc
        Nc = np.full(n, 2.81E25)

        global Nv
        Nv = np.full(n, 1.83E25)

        global Epsilon
        Epsilon = np.full(n, Epsilon_r * Epsilon_0)

        global mu_p
        mu_p = np.full(n, 0.045)

        global mu_n
        mu_n = np.full(n, 0.14)

        # Doping-Profile
        global C
        C = np.zeros(n)

        global Cau
        Cau = np.full(n, 0) # 1E-28
        
        global generation
        generation = np.full(n, 0.0)
 
        #
        global u
        u = np.zeros(3*n)

        #
        global u_old
        u_old = np.zeros(3*n)
        
        # 
        global b
        b = np.zeros(3*n)

        # 
        global A
        A = sparse.lil_matrix((3*n, 3*n))

        # Vector dx to be solved
        global x
        x = np.zeros(3*n)

        

def init_default_doping():
        """
        Definition of the default doping level.
        Left part doped with  :math:`N_\mathrm{a} = 1 \cdot 10^{24}\, \mathrm{m^{-3}}` and
        right part doped with :math:`N_\mathrm{d} = 1 \cdot 10^{24}\, \mathrm{m^{-3}}`
        
        """
        Na = 1E24
        Nd = 1E24

        for i in range(0,n):
                if i < n/2:
                    C[i] = -Na
                else:
                    C[i] = +Nd


        return None



init_geometry()
init_parameters()
init_default_doping()

