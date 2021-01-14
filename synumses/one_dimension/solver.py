import scipy.sparse as sparse 
from scipy.sparse.linalg import spsolve, qmr

import numpy as np

from synumses.one_dimension.scharfetter_gummel_bernoulli import bernoulli, jacobian, update_b, first_jacobian, first_update_b, transient_jacobian, transient_update_b

import synumses.one_dimension.parameters as parameters

from synumses.one_dimension.functions import calc_recombination, calc_p_density, calc_n_density, ohm_potential


def find_nearest(a, a0):
    "Index in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return idx

def time_step():
    
    parameters.u_old =  parameters.u


def solve_from_doping():
    
    for i in range(0,parameters.n):
            parameters.u[3*i + 0] = ohm_potential(parameters.C[i], parameters.Chi[i], parameters.Eg[i], parameters.Nc[i], parameters.Nv[i])
            parameters.u[3*i + 1] = 0.
            parameters.u[3*i + 2] = 0.
    return None



def solve_no_bias(norm_max = 1E-14, info = 'no'):
    
    solve_from_doping()
    
    while True:

        
        first_jacobian(0.0, 0.0)
        first_update_b(0.0, 0.0)

        calc_recombination()
        

        try:
            parameters.x = spsolve(parameters.A, parameters.b, use_umfpack=True)
        except RuntimeError:
            return False
        
        parameters.u = parameters.u + parameters.x
    
        norm = np.linalg.norm(parameters.x)

        if (info != 'no'):
            print("Norm of b: ", np.linalg.norm(parameters.b))
            print("Norm of x:" , np.linalg.norm(parameters.x))
        

        if (norm < norm_max):
            break

        
        #calc_p_density()
        #calc_n_density()
       
        
    return True


def solve_bias(Ua, Ub, norm_max  = 1E-6, info = 'no'):

    counter = 0

    while True:
        
        jacobian(Ua, Ub) 
        update_b(Ua, Ub)

        try:
            parameters.x = spsolve(parameters.A, parameters.b, use_umfpack=True)
        except RuntimeError:
            return False

        #parameters.x = qmr(parameters.A, parameters.b)
        parameters.u = parameters.u + parameters.x
    
        norm = np.linalg.norm(parameters.x)

        if (info != 'no'):
            print("Norm of b: ", np.linalg.norm(parameters.b))
            print("Norm of x:" , np.linalg.norm(parameters.x))

        #break
        

        if (np.isnan(norm) or (counter==20)):
            print("norm is nan")
            return False

        if (norm < norm_max):
            break

        counter = counter +1
  
    return True


def solve_bias_transient(Ua, Ub, dt, norm_max  = 1E-6, info = 'no'):
    
    parameters.dt = dt
    
    counter = 0

    while True: 
        
        transient_jacobian(Ua, Ub) 
        transient_update_b(Ua, Ub)

        try:
            parameters.x = spsolve(parameters.A, parameters.b, use_umfpack=True)
        except RuntimeError:
            return False

        #parameters.x = qmr(parameters.A, parameters.b)
        parameters.u = parameters.u + parameters.x
    
        norm = np.linalg.norm(parameters.x)

        if (info != 'no'):
            print("Norm of b: ", np.linalg.norm(parameters.b))
            print("Norm of x:" , np.linalg.norm(parameters.x))

        #break
        

        if (np.isnan(norm) or (counter==20)):
            print("norm is nan")
            return False

        if (norm < norm_max):
            break

        counter = counter +1


  
    return True

    
def solve_bias_center_boundary(Ua, Ub, Uc, Pos, doping, norm_max  = 1E-10, info = 'no'):
    
    prev_norm = 1E30

    while True:
    
        #calc_recombination()     
        
        jacobian(Ua, Ub) 
        update_b(Ua, Ub)
    
        base_pos =  find_nearest(parameters.pos_x, Pos)

        #
        # Boundary condition for potential
        #
        
        parameters.b[3*(base_pos+0)+0] =  ohm_potential(parameters.C[base_pos],
                                                        parameters.Chi[base_pos],
                                                        parameters.Eg[base_pos],
                                                        parameters.Nc[base_pos],
                                                        parameters.Nv[base_pos]) + Uc - parameters.u[3*(base_pos+0)+0] 

        parameters.A[3*(base_pos+0)+0, 3*(base_pos-1)+0] = 0
        parameters.A[3*(base_pos+0)+0, 3*(base_pos-1)+1] = 0
        parameters.A[3*(base_pos+0)+0, 3*(base_pos-1)+2] = 0

        parameters.A[3*(base_pos+0)+0, 3*(base_pos+0)+0] = 1
        parameters.A[3*(base_pos+0)+0, 3*(base_pos+0)+1] = 0
        parameters.A[3*(base_pos+0)+0, 3*(base_pos+0)+2] = 0

        parameters.A[3*(base_pos+0)+0, 3*(base_pos+1)+0] = 0
        parameters.A[3*(base_pos+0)+0, 3*(base_pos+1)+1] = 0
        parameters.A[3*(base_pos+0)+0, 3*(base_pos+1)+2] = 0

        if (doping =="p"):
            #
            # Boundary condition for hole fermi level
            #
            parameters.b[3*(base_pos+0)+1] =  Uc - parameters.u[3*(base_pos+0)+1] 

            parameters.A[3*(base_pos+0)+1, 3*(base_pos-1)+0] = 0
            parameters.A[3*(base_pos+0)+1, 3*(base_pos-1)+1] = 0
            parameters.A[3*(base_pos+0)+1, 3*(base_pos-1)+2] = 0

            parameters.A[3*(base_pos+0)+1, 3*(base_pos+0)+0] = 0
            parameters.A[3*(base_pos+0)+1, 3*(base_pos+0)+1] = 1
            parameters.A[3*(base_pos+0)+1, 3*(base_pos+0)+2] = 0

            parameters.A[3*(base_pos+0)+1, 3*(base_pos+1)+0] = 0
            parameters.A[3*(base_pos+0)+1, 3*(base_pos+1)+1] = 0
            parameters.A[3*(base_pos+0)+1, 3*(base_pos+1)+2] = 0

        elif (doping =="n"):
            #
            # Boundary condition for electron fermi level
            #
            parameters.b[3*(base_pos+0)+2] =  Uc - parameters.u[3*(base_pos+0)+2] 

            parameters.A[3*(base_pos+0)+2, 3*(base_pos-1)+0] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos-1)+1] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos-1)+2] = 0

            parameters.A[3*(base_pos+0)+2, 3*(base_pos+0)+0] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos+0)+1] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos+0)+2] = 1

            parameters.A[3*(base_pos+0)+2, 3*(base_pos+1)+0] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos+1)+1] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos+1)+2] = 0

        try:
            parameters.x = spsolve(parameters.A, parameters.b, use_umfpack=True)
        except RuntimeError:
            return False

        #parameters.x = qmr(parameters.A, parameters.b)
        parameters.u = parameters.u + parameters.x
    
        norm = np.linalg.norm(parameters.x)

        if (info != 'no'):
            print("Norm of b: ", np.linalg.norm(parameters.b))
            print("Norm of x:" , np.linalg.norm(parameters.x))
        
        if (np.isnan(norm)):
            print("norm is nan")
            return False

        if ((prev_norm <= norm) and (norm < norm_max)):
            break
             
        prev_norm = norm
  
    return True


def solve_bias_center_boundary_transient(Ua, Ub, Uc, Pos, doping, dt, norm_max  = 1E-10, info = 'no'):

    prev_norm = 1E30

    parameters.dt = dt
    
    while True:
    
        #calc_recombination()     
        
        transient_jacobian(Ua, Ub) 
        transient_update_b(Ua, Ub)
    
        base_pos =  find_nearest(parameters.pos_x, Pos)

        #
        # Boundary condition for potential
        #
        
        parameters.b[3*(base_pos+0)+0] =  ohm_potential(parameters.C[base_pos],
                                                        parameters.Chi[base_pos],
                                                        parameters.Eg[base_pos],
                                                        parameters.Nc[base_pos],
                                                        parameters.Nv[base_pos]) + Uc - parameters.u[3*(base_pos+0)+0] 

        parameters.A[3*(base_pos+0)+0, 3*(base_pos-1)+0] = 0
        parameters.A[3*(base_pos+0)+0, 3*(base_pos-1)+1] = 0
        parameters.A[3*(base_pos+0)+0, 3*(base_pos-1)+2] = 0

        parameters.A[3*(base_pos+0)+0, 3*(base_pos+0)+0] = 1
        parameters.A[3*(base_pos+0)+0, 3*(base_pos+0)+1] = 0
        parameters.A[3*(base_pos+0)+0, 3*(base_pos+0)+2] = 0

        parameters.A[3*(base_pos+0)+0, 3*(base_pos+1)+0] = 0
        parameters.A[3*(base_pos+0)+0, 3*(base_pos+1)+1] = 0
        parameters.A[3*(base_pos+0)+0, 3*(base_pos+1)+2] = 0

        if (doping =="p"):
            #
            # Boundary condition for hole fermi level
            #
            parameters.b[3*(base_pos+0)+1] =  Uc - parameters.u[3*(base_pos+0)+1] 

            parameters.A[3*(base_pos+0)+1, 3*(base_pos-1)+0] = 0
            parameters.A[3*(base_pos+0)+1, 3*(base_pos-1)+1] = 0
            parameters.A[3*(base_pos+0)+1, 3*(base_pos-1)+2] = 0

            parameters.A[3*(base_pos+0)+1, 3*(base_pos+0)+0] = 0
            parameters.A[3*(base_pos+0)+1, 3*(base_pos+0)+1] = 1
            parameters.A[3*(base_pos+0)+1, 3*(base_pos+0)+2] = 0

            parameters.A[3*(base_pos+0)+1, 3*(base_pos+1)+0] = 0
            parameters.A[3*(base_pos+0)+1, 3*(base_pos+1)+1] = 0
            parameters.A[3*(base_pos+0)+1, 3*(base_pos+1)+2] = 0

        elif (doping =="n"):
            #
            # Boundary condition for electron fermi level
            #
            parameters.b[3*(base_pos+0)+2] =  Uc - parameters.u[3*(base_pos+0)+2] 

            parameters.A[3*(base_pos+0)+2, 3*(base_pos-1)+0] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos-1)+1] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos-1)+2] = 0

            parameters.A[3*(base_pos+0)+2, 3*(base_pos+0)+0] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos+0)+1] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos+0)+2] = 1

            parameters.A[3*(base_pos+0)+2, 3*(base_pos+1)+0] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos+1)+1] = 0
            parameters.A[3*(base_pos+0)+2, 3*(base_pos+1)+2] = 0

        try:
            parameters.x = spsolve(parameters.A, parameters.b, use_umfpack=True)
        except RuntimeError:
            return False

        #parameters.x = qmr(parameters.A, parameters.b)
        parameters.u = parameters.u + parameters.x
    
        norm = np.linalg.norm(parameters.x)

        if (info != 'no'):
            print("Norm of b: ", np.linalg.norm(parameters.b))
            print("Norm of x:" , np.linalg.norm(parameters.x))
        
        if (np.isnan(norm)):
            print("norm is nan")
            return False

        if ((prev_norm <= norm) and (norm < norm_max)):
            break
             
        prev_norm = norm
  
    return True


