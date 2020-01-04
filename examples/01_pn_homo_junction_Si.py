#!/usr/bin/env python
# coding: utf-8

# # Simulation of a Si homo pn-junction

# The following code loads the needed modules.

# In[1]:


import numpy as np

import matplotlib.pyplot as plt

from synumses.one_dimension.scharfetter_gummel_bernoulli import bernoulli, jacobian, update_b, first_jacobian, first_update_b, hole_current_density, electron_current_density

import synumses.one_dimension.parameters as parameters

from synumses.one_dimension.functions import calc_p_density, calc_n_density


import synumses.one_dimension.solver as solver


# # Simulation of the default device
# 
# By default the package uses the following parameters for silicon:<br>
# * **parameters.n  = 400**, 400 cell points
# * **parameters.lx = 400E-9**, total length of device $l_\mathrm{x} = 400\,\mathrm{nm}$
# * **parameters.C[000:199]   = -1E24**, left part of device is p-doped $N_\mathrm{A} = 10^{24}\,\mathrm{m}^{-3}$
# * **parameters.C[200:399] = +1E24**, right part of device is p-doped $N_\mathrm{D} = 10^{24}\,\mathrm{m}^{-3}$
# * **parameters.Epsilon[0:399] = $\varepsilon_0 \cdot \varepsilon_\mathrm{r, Si}$**, materialparameter von silicon 
# * **parameters.Ec[0:399] = $1.12\,\mathrm{eV}$**, energie level of conduction band
# * **parameters.Ev[0:399] = $0\,\mathrm{eV}$**, energie level of valence band
# * **parameters.Cau[0:399] = $2.3\cdot 10^{-20} \mathrm{\dfrac{m^3}{s}}$** defines the recombination rate :
#     $\dfrac{\mathrm{d}n}{\mathrm{d}t} = Cau \left(n \cdot p - n_\mathrm{i}^2 \right)$
# 

# ## Plot some default parameters
# * Band diagram
# * Doping

# In[2]:


#
# Plot electron affinity  and \n band gap
#
fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Electron affinity $\chi$ and \n band gap", y= 1.05, fontsize=16)

axis.plot(parameters.pos_x * 1E9, parameters.Chi, label='Electron affinity $\chi$')
axis.plot(parameters.pos_x * 1E9, parameters.Eg, label='Band gap')


plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"Energy in $\mathrm{eV}$")


axis.legend(loc = 6)
plt.show()

#
# Plot band diagram
#
fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Electron affinity $\chi$ and \n band gap", y= 1.05, fontsize=16)

axis.plot(parameters.pos_x * 1E9, -parameters.Chi, label='Conduction band')
axis.plot(parameters.pos_x * 1E9, -parameters.Eg - parameters.Chi, label='Valence band')


plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"Energy in $\mathrm{eV}$")


axis.legend(loc = 6)
plt.show()

#
# Plot doping profile of the pn-junction
#
fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Doping profile of the pn-junction", fontsize=16)

axis.plot(parameters.pos_x * 1E9, parameters.C, label='Doping')


plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"Doping in $\mathrm{\dfrac{1}{m^3}}$")
plt.show()


# ## Calculate the potential according the doping
# The following function calculates the potential according the doping level.
# This is needed the have a first guess of the potential.

# In[3]:


solver.solve_from_doping()


# Now, we plot the potential according the doping level. It's not smooth.

# In[4]:


fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Potential from Doping", fontsize=16)

plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"Energy in $\mathrm{eV}$")

axis.plot(parameters.pos_x * 1E9, parameters.u[0::3])
plt.show()


# ## Calculate the potential
# Now, let's calculate the potential considering diffusion but no biasing.<br>
# One simulation using **solver.solve_no_bias()** must be performed befor simulatoins with biasing.

# In[5]:


solver.solve_no_bias()


# ## Plot the potential

# In[6]:


fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Potential", fontsize=16)

plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"Energy in $\mathrm{eV}$")

axis.plot(parameters.pos_x * 1E9, parameters.u[0::3])
plt.show()


# ## Calculate the potential and quasi Fermi levels without biasing

# In[7]:


solver.solve_bias(0,0)


# ## Plot the results
# * potential,
# * electron and hole density,
# * electron and hole current density, and

# In[8]:


#
# Plot potential
#
fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Potential", fontsize=16)

plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"$\varphi$ in $\mathrm{eV}$")

axis.plot(parameters.pos_x * 1E9, parameters.u[0::3], label='Potential')
axis.legend()
plt.show()

#
# Plot band diagram
#
fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Potential", fontsize=16)

plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"$\varphi$ in $\mathrm{eV}$")

axis.plot(parameters.pos_x * 1E9, -parameters.u[0::3] - parameters.Chi -  parameters.Eg, label='Valence band')
axis.plot(parameters.pos_x * 1E9, -parameters.u[0::3] - parameters.Chi                 , label='Conduction band')
axis.plot(parameters.pos_x * 1E9, -parameters.u[1::3]                                  , label='Quasi-Fermi-level of holes')
axis.plot(parameters.pos_x * 1E9, -parameters.u[2::3]                                  , label='Quasi-Fermi-level of electrons')

axis.legend()
plt.show()


#
# Plot electron and hole density
#
p = calc_p_density()
n = calc_n_density()

fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Density of holes and electrons", fontsize=16)

plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"n in $\mathrm{m^{-3}}$")

axis.semilogy(parameters.pos_x * 1E9, p, label='Hole density')
axis.semilogy(parameters.pos_x * 1E9, n, label='Electron density')
axis.legend()
plt.show()

#
# Plot space charge
#
p = calc_p_density()
n = calc_n_density()

fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Space charge", fontsize=16)

plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"n in $\mathrm{C \cdot m^{-3}}$")

axis.plot(parameters.pos_x * 1E9, p-n+parameters.C, label='Space charge')

axis.legend()
plt.show()

#
# Plot electron and hole current density 
#
j_p =     hole_current_density()
j_n = electron_current_density()

fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Current density of holes and electrons", fontsize=16)

plt.xlabel(r"$x$ in $\mathrm{nm}$")
plt.ylabel(r"$j$ in $\mathrm{\dfrac{A}{m^{-2}}}$")

axis.plot(parameters.pos_x * 1E9, j_p, label='Hole current density')
axis.plot(parameters.pos_x * 1E9, j_n, label='Electron current density')
axis.legend()
plt.show()


# ## Simulating a voltage sweep
# Now let's simulate a voltage sweep to get the diode characteristic. <br>
# The bias points are printed and stored in the arrays **voltage** and **current_density**.

# In[9]:


u_start = 0.0
u_stop  = 0.7
u_step  = 0.025

voltage = []
current_density = []

bias_points = np.linspace(u_start, u_stop, int((u_stop-u_start)/(u_step))+2)
print(bias_points)
for bias_point in bias_points:

    solver.solve_bias(bias_point,0)   
    
    j = np.mean(hole_current_density() + electron_current_density())
    voltage.append(bias_point)
    current_density.append(j)
    print(bias_point, ",", j)
  


# In[10]:


#
# Diode charactersistics
#
fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Diode characteristics", fontsize=16)

axis.plot(voltage, current_density, label='Current density')

axis.set_xlabel(r"Voltage in $\mathrm{V}$")
axis.set_ylabel(r"Current density in $\mathrm{\dfrac{A}{m^2}}$")

axis.legend()
#plt.ylim(-1E-6,1E-6)
plt.show()

#
# Plot potential and quasi-Fermi-levels
#
fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Conduction band, valence band and \nquasi Fermi levels", y= 1.05, fontsize=16)

plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"Energy in $\mathrm{eV}$")

axis.plot(parameters.pos_x * 1E9, -parameters.u[0::3] - parameters.Chi -  parameters.Eg, label='Valence band')
axis.plot(parameters.pos_x * 1E9, -parameters.u[0::3] - parameters.Chi                 , label='Conduction band')
axis.plot(parameters.pos_x * 1E9, -parameters.u[1::3]                                  , label='Quasi-Fermi-level of holes')
axis.plot(parameters.pos_x * 1E9, -parameters.u[2::3]                                  , label='Quasi-Fermi-level of electrons')

plt.legend()
plt.show()

#
# Plot electron and hole current density 
#
j_p =     hole_current_density()
j_n = electron_current_density()

fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Current density of holes and electrons", fontsize=16)

plt.xlabel(r"$x$ in $\mathrm{nm}$")
plt.ylabel(r"$j$ in $\mathrm{\dfrac{A}{m^{-2}}}$")

axis.plot(parameters.pos_x * 1E9, j_p, label='Hole current density')
axis.plot(parameters.pos_x * 1E9, j_n, label='Electron current density')
axis.legend()
plt.show()


# Electron and hole density outside the space cahreg

p = calc_p_density()
n = calc_n_density()

fig, axis = plt.subplots(1,1, sharey=True)
fig.suptitle("Density of holes and electrons \n outside the space charge \n without recombination",  y=1.10, fontsize=16)

plt.xlabel(r"x in $\mathrm{nm}$")
plt.ylabel(r"p/n in $\mathrm{m^{-3}}$")

#axis.semilogy(p, label='Hole density')
axis.plot(parameters.pos_x * 1E9, p, label='Hole density')
#axis.semilogy(n, label='Electron density')
axis.plot(parameters.pos_x * 1E9, n, label='Electron density')

axis.legend()
plt.ylim(0, 1E20)
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




