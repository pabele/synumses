import scipy.sparse as sparse 
from scipy.sparse.linalg import spsolve

import numpy as np

from modells.scharfetter_gummel_bernoulli import bernoulli, jacobian, update_b, first_jacobian, first_update_b

import modells.parameters as parameters

from modells.functions import calc_recombination

def solve_no_bias():
    
    parameters.init_potential()
    
    while True:

        first_jacobian(0.0, 0.0)
        first_update_b(0.0, 0.0)

        
        parameters.x = spsolve(parameters.A, parameters.b)
        parameters.u = parameters.u + parameters.x
    
        norm = np.linalg.norm(parameters.b)
        #print("Norm of b:")
        #print(norm)
        
        if (norm < 1.0E-14):
            break

    return None

def solve_bias(Ua, Ub):
    
    prev_norm = 1E300
    
    while True:
    
        calc_recombination()     
        
        jacobian(Ua, Ub) 
        update_b(Ua, Ub)
    
        parameters.x = spsolve(parameters.A, parameters.b)        
        parameters.u = parameters.u + parameters.x
    
        norm = np.linalg.norm(parameters.b)
        #print("Norm of b:")
        #print(norm)
        
        if (prev_norm < norm):
            break
             
        prev_norm = norm
    return None
