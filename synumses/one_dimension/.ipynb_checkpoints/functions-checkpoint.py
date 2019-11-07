import numpy as np

import modells.parameters as parameters


def calc_p_density():
    
    parameters.p_density = parameters.Nv*np.exp(parameters.q*( parameters.u[1::3] - parameters.u[0::3] + parameters.Ev)/(parameters.kB*parameters.T))
          
    return None



def hole_current_density():
    from modells.scharfetter_gummel_bernoulli import bernoulli
    j_p = np.zeros(parameters.n)

    for i in range(0,parameters.n-1):
       
        j_p[i] =(
            +(parameters.q*parameters.mu_p[i]*parameters.Ut)/(parameters.dx)*(
            +bernoulli(+(parameters.u[3*(i+1)+0]-parameters.u[3*(i+0)+0])/(parameters.Ut)) * parameters.Nv[i]*np.exp(parameters.q*( parameters.u[3*(i+0)+1] - parameters.u[3*(i+0)+0] + (parameters.Ev[i] + parameters.Ev[i+1])/2.)/(parameters.kB*parameters.T))
            -bernoulli(-(parameters.u[3*(i+1)+0]-parameters.u[3*(i+0)+0])/(parameters.Ut)) * parameters.Nv[i+1]*np.exp(parameters.q*( parameters.u[3*(i+1)+1] - parameters.u[3*(i+1)+0] + (parameters.Ev[i] + parameters.Ev[i+1])/2.)/(parameters.kB*parameters.T)))
        )

    i = parameters.n-1
    j_p[i] =  j_p[i-1]
#    j_p[i] = (
#        +(parameters.q*parameters.mu_p[i]*parameters.Ut)/(parameters.dx)*(
#        +bernoulli((+(ohm_potential(parameters.C[parameters.n-1], parameters.Ec[parameters.n-1], parameters.Ev[parameters.n-1], parameters.Nc[parameters.n-1], parameters.Nv[parameters.n-1]) + Ub) -parameters.u[3*(i+0)+0])/(parameters.Ut)) * parameters.Nv[i]*np.exp(parameters.q*( parameters.u[3*(i+0)+1] - parameters.u[3*(i+0)+0] + parameters.Ev[i])/(parameters.kB*parameters.T))
#        -bernoulli((-(ohm_potential(parameters.C[parameters.n-1], parameters.Ec[parameters.n-1], parameters.Ev[parameters.n-1], parameters.Nc[parameters.n-1], parameters.Nv[parameters.n-1]) + Ub) -parameters.u[3*(i+0)+0])/(parameters.Ut)) * parameters.Nv[i]*np.exp(parameters.q*( Ub - ( ohm_potential(parameters.C[parameters.n-1], parameters.Ec[parameters.n-1], parameters.Ev[parameters.n-1], parameters.Nc[parameters.n-1], parameters.Nv[parameters.n-1]) + Ub)   + parameters.Ev[i])/(parameters.kB*parameters.T)))
#    )

    
    return j_p


def calc_n_density():

    parameters.n_density = parameters.Nc*np.exp(parameters.q*( parameters.u[0::3] - parameters.u[2::3] - parameters.Ec)/(parameters.kB*parameters.T))
   
    return None



def electron_current_density():
    from modells.scharfetter_gummel_bernoulli import bernoulli
    j_n = np.zeros( parameters.n)

    
    for i in range(0, parameters.n-1):
       
        j_n[i] =(
            -( parameters.q* parameters.mu_n[i]* parameters.Ut)/( parameters.dx)*(
            +bernoulli(-( parameters.u[3*(i+1)+0]- parameters.u[3*(i+0)+0])/( parameters.Ut)) *  parameters.Nc[i]*np.exp( parameters.q*(  parameters.u[3*(i+0)+0] -  parameters.u[3*(i+0)+2] -  (parameters.Ec[i] + parameters.Ec[i+1])/2.)/( parameters.kB* parameters.T))
            -bernoulli(+( parameters.u[3*(i+1)+0]- parameters.u[3*(i+0)+0])/( parameters.Ut)) *  parameters.Nc[i+1]*np.exp( parameters.q*(  parameters.u[3*(i+1)+0] -  parameters.u[3*(i+1)+2] -  (parameters.Ec[i] + parameters.Ec[i+1])/2.)/( parameters.kB* parameters.T)))
        )

    i =  parameters.n-1
    j_n[i] =  j_n[i-1]
#    j_n[i] = (
#        -( parameters.q*parameters.mu_n[i]* parameters.Ut)/( parameters.dx)*(
#        +bernoulli(-( ohm_potential(parameters.C[i], parameters.Ec[i], parameters.Ev[i], parameters.Nc[i], parameters.Nv[i]) + Ub) -parameters.u[3*(i+0)+0])/(parameters.Ut) * parameters.Nc[i]*np.exp(parameters.q*( parameters.u[3*(i+0)+0] - parameters.u[3*(i+0)+2] - parameters.Ec[i])/(parameters.kB*parameters.T))
#        -bernoulli(+(ohm_potential(parameters.C[i], parameters.Ec[i], parameters.Ev[i], parameters.Nc[i], parameters.Nv[i]) + Ub) -parameters.u[3*(i+0)+0])/(parameters.Ut) * parameters.Nc[i]*np.exp(parameters.q*( (ohm_potential(parameters.C[i], parameters.Ec[i], parameters.Ev[i], parameters.Nc[i], parameters.Nv[i]) + Ub) -Ub    - parameters.Ec[i])/(parameters.kB*parameters.T)))
#    )

    
    return j_n



def ohm_potential(C, Ec, Ev, Nc, Nv):

    from modells.parameters import q, kB, Ut, T
        
    ni = np.sqrt(Nc*Nv)*np.exp(-((Ec-Ev)*q)/(2.*kB*T))
    Phi = (Ec + Ev)/2. - 0.5*Ut*np.log(Nc/Nv) + Ut*np.arcsinh(C/(2.*ni))
    #print(Phi)
    
    #if (C > 0):
    #    Phi = Ut*np.log((np.sqrt(C**2 + 4.*ni**2)+C)/(2.*ni))
    #else:
    #    Phi = -Ut*np.log((np.sqrt(C**2 + 4.*ni**2)-C)/(2.*ni))
    # 
    #print(Phi)
    
    return(Phi)


def calc_ni_density():

    parameters.ni_density = (
        np.sqrt(parameters.Nc*parameters.Nv)
        *
        np.exp(-((parameters.Ec-parameters.Ev)*parameters.q)
               /
               (2.*parameters.kB*parameters.T))
    )

    return None

def calc_recombination():
    
    calc_p_density()
    calc_n_density()
    #calc_ni_density()
    
    parameters.recombination = (
        parameters.Cau
        *
        (parameters.p_density * parameters.n_density -
         parameters.ni_density**2
        )
    )
    
#    parameters.recombination = (
#        (parameters.p_density * parameters.n_density -
#         parameters.ni_density**2
#        )
#        /
#        (parameters.tau_p*parameters.p_density + 
#         parameters.tau_n*parameters.n_density)
#    )

    return None

def calc_generation(I_0, alpha):
    
    return None
    
