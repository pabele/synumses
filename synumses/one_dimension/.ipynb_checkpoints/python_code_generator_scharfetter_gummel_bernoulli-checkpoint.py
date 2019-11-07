from sympy import symbols, Symbol 
from sympy import Function
from sympy import solve, pretty
from sympy import diff
from sympy import sin, cos, exp, ln, sinh, cosh
from sympy import simplify, trigsimp

import time

def bernoulli(expr):
    #return 1.-0.5*expr
    return 1.0 - expr/2.0 + expr**2/12.0 - expr**4/720.0 + expr**6/30240.0
    #return expr/(exp(expr) - 1)


Phi_n_k, Phi_n_l = symbols('Phi_n_k, Phi_n_l')
Phi_p_k, Phi_p_l = symbols('Phi_p_k, Phi_p_l')
Psi_k, Psi_l     = symbols('Psi_k, Psi_l')

Ev_k, Ev_l = symbols('Ev_k, Ev_l')
Nv_k, Nv_l = symbols('Nv_k, Nv_l')
Ec_k, Ec_l = symbols('Ec_k, Ec_l')
Nc_k, Nc_l = symbols('Nc_k, Nc_l')

Psi_m1, Psi_00, Psi_p1       = symbols('parameters.u[3*(i-1)+0], parameters.u[(3*i)+0], parameters.u[3*(i+1)+0]')
Phi_p_m1, Phi_p_00, Phi_p_p1 = symbols('parameters.u[3*(i-1)+1], parameters.u[(3*i)+1], parameters.u[3*(i+1)+1]')
Phi_n_m1, Phi_n_00, Phi_n_p1 = symbols('parameters.u[3*(i-1)+2], parameters.u[(3*i)+2], parameters.u[3*(i+1)+2]')


Ev_m1, Ev_00, Ev_p1       = symbols('parameters.Ev[i-1], parameters.Ev[i+0], parameters.Ev[i+1]')
Nv_m1, Nv_00, Nv_p1       = symbols('parameters.Nv[i-1], parameters.Nv[i+0], parameters.Nv[i+1]')
Ec_m1, Ec_00, Ec_p1       = symbols('parameters.Ec[i-1], parameters.Ec[i+0], parameters.Ec[i+1]')
Nc_m1, Nc_00, Nc_p1       = symbols('parameters.Nc[i-1], parameters.Nc[i+0], parameters.Nc[i+1]')

dx = symbols('parameters.dx')




C_00 = symbols('parameters.C[i]')
Epsilon_00 = symbols('parameters.Epsilon[i]')

np_exp = Function('np.exp')

kB = symbols('parameters.kB')
T  = symbols('parameters.T')
q  = symbols('parameters.q')

Ut = (kB * T) / q

mu_n, mu_p = symbols('parameters.mu_n[i], parameters.mu_p[i]')

p = Function('p')
n = Function('n')

p = Nv_00*exp(q*( Ev_00 + Phi_p_00 - Psi_00)/(kB*T))
n = Nc_00*exp(q*(-Ec_00 - Phi_n_00 + Psi_00)/(kB*T))


#poisson = ((Psi_p1 - 2*Psi_00 + Psi_m1) + q / Epsilon * (C_00 + p -n)*dx**2).subs([(Psi_k, Psi_00),(Psi_l, Psi_p1),(Phi_p_k, Phi_p_00),(Phi_p_l,Phi_p_p1)])
poisson = ((Psi_p1 - 2*Psi_00 + Psi_m1) + q / Epsilon_00 * (C_00 + p -n)*dx**2)

j_p =+(q*mu_p*Ut)/(dx)*(
    +bernoulli(+(Psi_l-Psi_k)/(Ut)) * Nv_k*exp(q*(Phi_p_k - Psi_k + Ev_k)/(kB*T))
    -bernoulli(-(Psi_l-Psi_k)/(Ut)) * Nv_l*exp(q*(Phi_p_l - Psi_l + Ev_l)/(kB*T))
)

j_n = -(q*mu_n*Ut)/(dx)*(
    +bernoulli(-(Psi_l-Psi_k)/(Ut)) * Nc_k*exp(q*(Psi_k - Phi_n_k - Ec_k)/(kB*T))
    -bernoulli(+(Psi_l-Psi_k)/(Ut)) * Nc_l*exp(q*(Psi_l - Phi_n_l - Ec_l)/(kB*T))
)

div_j_p = (
    + j_p.subs([(Psi_k, Psi_00),(Psi_l, Psi_p1),(Phi_p_k, Phi_p_00),(Phi_p_l,Phi_p_p1),(Ev_k, Ev_00),(Ev_l, Ev_p1),(Nv_k, Nv_00),(Nv_l, Nv_p1)])
    - j_p.subs([(Psi_k, Psi_m1),(Psi_l, Psi_00),(Phi_p_k, Phi_p_m1),(Phi_p_l,Phi_p_00),(Ev_k, Ev_m1),(Ev_l, Ev_00),(Nv_k, Nv_m1),(Nv_l, Nv_00)])
)

div_j_n = (
    + j_n.subs([(Psi_k, Psi_00),(Psi_l, Psi_p1),(Phi_n_k, Phi_n_00),(Phi_n_l,Phi_n_p1),(Ec_k, Ec_00),(Ec_l, Ec_p1),(Nc_k, Nc_00),(Nc_l, Nc_p1)])
    - j_n.subs([(Psi_k, Psi_m1),(Psi_l, Psi_00),(Phi_n_k, Phi_n_m1),(Phi_n_l,Phi_n_00),(Ec_k, Ec_m1),(Ec_l, Ec_00),(Nc_k, Nc_m1),(Nc_l, Nc_00)])
)



#############################
# Start of module generation
#############################



functions = [["3*i+0",poisson, "poisson"],
             ["3*i+1",div_j_p, "div_j_p"],
             ["3*i+2",div_j_n, "div_j_n"]]


partial_derivatives = [["3*(i-1)+0",Psi_m1  ],
                      ["3*(i+0)+0",Psi_00  ],
                      ["3*(i+1)+0",Psi_p1  ],
                      ["3*(i-1)+1",Phi_p_m1],
                      ["3*(i+0)+1",Phi_p_00],
                      ["3*(i+1)+1",Phi_p_p1],
                      ["3*(i-1)+2",Phi_n_m1],
                      ["3*(i+0)+2",Phi_n_00],
                      ["3*(i+1)+2",Phi_n_p1]
]



substitudes = {"left"  : ((Psi_m1, '(U_Ohm(parameters.C[0             ], parameters.Ec[0             ], parameters.Ev[0             ], parameters.Nc[0             ], parameters.Nv[0             ]) + Ua)'),
                          (Phi_p_m1, 'Ua'),
                          (Phi_n_m1, 'Ua'),
                          (exp,np_exp)),
               "right" : ((Psi_p1, '(U_Ohm(parameters.C[parameters.n-1], parameters.Ec[parameters.n-1], parameters.Ev[parameters.n-1], parameters.Nc[parameters.n-1], parameters.Nv[parameters.n-1]) + Ub)'),
                          (Phi_p_p1, 'Ub'),
                          (Phi_n_p1, 'Ub'),
                          (Ev_p1, Ev_00),
                          (Nv_p1, Nv_00),
                          (Ec_p1, Ec_00),
                          (Nc_p1, Nc_00),
                          (exp,np_exp)),
               "center": ((1,1),(exp,np_exp))
    }

print("####################################")
print("### This code was automatically  ###")
print("###        generated             ###")
print("###    on",time.strftime("%d.%m.%Y at %H:%M"), "   ###")
print("####################################")

print()

print("#################################")
print("###########  Import   ###########")
print("#################################")
print("import numpy as np")
print()
print("import modells.parameters as parameters")
print("from modells.contacts import U_Ohm")

print()
print()

print("#################################")
print("########### Bernoulli ###########")
print("#################################")
print("def bernoulli(x):")

x = Symbol('x')
print("\t return",bernoulli(x))
print()
print()


print("#################################")
print("########### Function ###########")
print("#################################")

print("def update_b(Ua, Ub):")
print()
print("\t for i in range(0, parameters.n):")

for s in ["left", "right", "center"]:
    if (s =="left"):
        print("\t \t if (i==0) :")
    elif (s=="right"):
        print("\t \t elif(i==parameters.n-1):")
    else:
        print("\t \t else:")

    print("\t \t \t #################")
    print("\t \t \t ### ",s, "###")
    print("\t \t \t #################")
    for function in functions:
        
        print("\t \t \t parameters.b["+function[0]+"] = ", function[1].subs(substitudes[s]))

print() 
print("\t return None")
print()
print()
print()

print("#################################")
print("###########  Jacobi  ###########")
print("#################################")

print("def jacobian(Ua, Ub):")
print()
print("\t for i in range(0, parameters.n):")

for s in ["left", "right", "center"]:

    if (s =="left"):
        print("\t \t if (i==0) :")
    elif (s=="right"):
        print("\t \t elif(i==parameters.n-1):")
    else:
        print("\t \t else:")

    print()
    print("\t \t \t #################")
    print("\t \t \t ### ",s, "###")
    print("\t \t \t #################")

    #print("# ",substitudes[s])

    for function in functions:
        print
        print("\t \t \t #######################")
        print("\t \t \t ### ",function[2]," ###")
        print("\t \t \t #######################")
        for partial_derivative in partial_derivatives:
            if ((s == "left"  and (partial_derivative[1] in [Psi_m1, Phi_p_m1, Phi_n_m1])) or
                (s == "right" and (partial_derivative[1] in [Psi_p1, Phi_p_p1, Phi_n_p1]))):
                print("\t \t \t ###"+"parameters.A["+function[0]+","+partial_derivative[0]+"]", " = -(", diff(function[1], partial_derivative[1]).subs(substitudes[s]), ")")
            else:
                print("\t \t \t parameters.A["+function[0]+","+partial_derivative[0]+"]", " = -(", diff(function[1], partial_derivative[1]).subs(substitudes[s]), ")")
           
            
print()
print("\t return ")
print()
print()

poisson =  ((Psi_p1 - 2*Psi_00 + Psi_m1) + q / Epsilon_00 * (C_00 + p -n)*dx**2)
div_j_p = Phi_p_00
div_j_n = Phi_n_00

functions = [["3*i+0",poisson, "poisson"],
             ["3*i+1",div_j_p, "div_j_p"],
             ["3*i+2",div_j_n, "div_j_n"]]

partial_derivatives = [["3*(i-1)+0",Psi_m1  ],
                      ["3*(i+0)+0",Psi_00  ],
                      ["3*(i+1)+0",Psi_p1  ],
                      ["3*(i-1)+1",Phi_p_m1],
                      ["3*(i+0)+1",Phi_p_00],
                      ["3*(i+1)+1",Phi_p_p1],
                      ["3*(i-1)+2",Phi_n_m1],
                      ["3*(i+0)+2",Phi_n_00],
                      ["3*(i+1)+2",Phi_n_p1]
]

substitudes = {"left"  : ((Psi_m1, '(U_Ohm(parameters.C[0             ], parameters.Ec[0             ], parameters.Ev[0             ], parameters.Nc[0             ], parameters.Nv[0             ]) + Ua)'),
                          (Phi_p_m1, 'Ua'),
                          (Phi_n_m1, 'Ua'),
                          (exp,np_exp)),
               "right" : ((Psi_p1, '(U_Ohm(parameters.C[parameters.n-1], parameters.Ec[parameters.n-1], parameters.Ev[parameters.n-1], parameters.Nc[parameters.n-1], parameters.Nv[parameters.n-1]) + Ub)'),
                          (Phi_p_p1, 'Ub'),
                          (Phi_n_p1, 'Ub'),
                          (Ev_p1, Ev_00),
                          (Nv_p1, Nv_00),
                          (Ec_p1, Ec_00),
                          (Nc_p1, Nc_00),
                          (exp,np_exp)),
               "center": ((1,1),(exp,np_exp))
    }


print("#######################################")
print("########### First Function ###########")
print("#######################################")

print("def first_update_b(Ua, Ub):")
print()
print("\t for i in range(0, parameters.n):")

for s in ["left", "right", "center"]:
    if (s =="left"):
        print("\t \t if (i==0) :")
    elif (s=="right"):
        print("\t \t elif(i==parameters.n-1):")
    else:
        print("\t \t else:")

    print("\t \t \t #################")
    print("\t \t \t ### ",s, "###")
    print("\t \t \t #################")
    for function in functions:
        
        print("\t \t \t parameters.b["+function[0]+"] = ", function[1].subs(substitudes[s]))

print() 
print("\t return None")
print()
print()
print()

print("######################################")
print("###########  First Jacobi  ###########")
print("######################################")

print("def first_jacobian(Ua, Ub):")
print()
print("\t for i in range(0, parameters.n):")

for s in ["left", "right", "center"]:

    if (s =="left"):
        print("\t \t if (i==0) :")
    elif (s=="right"):
        print("\t \t elif(i==parameters.n-1):")
    else:
        print("\t \t else:")

    print()
    print("\t \t \t #################")
    print("\t \t \t ### ",s, "###")
    print("\t \t \t #################")

    #print("# ",substitudes[s])

    for function in functions:
        print
        print("\t \t \t #######################")
        print("\t \t \t ### ",function[2]," ###")
        print("\t \t \t #######################")
        for partial_derivative in partial_derivatives:
            if ((s == "left"  and (partial_derivative[1] in [Psi_m1, Phi_p_m1, Phi_n_m1])) or
                (s == "right" and (partial_derivative[1] in [Psi_p1, Phi_p_p1, Phi_n_p1]))):
                print("\t \t \t ###"+"parameters.A["+function[0]+","+partial_derivative[0]+"]", " = -(", diff(function[1], partial_derivative[1]).subs(substitudes[s]), ")")
            else:
                print("\t \t \t parameters.A["+function[0]+","+partial_derivative[0]+"]", " = -(", diff(function[1], partial_derivative[1]).subs(substitudes[s]), ")")
           
            
print()
print("\t return ")
print()
print()
