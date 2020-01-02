"""
This module generates the python code for the
module **scharfetter_gummel_bernoulli.py**.
It solves the *Van Roosbroeck equations* using the 
*Scharfetter Gummel Scheme*.

"""

from sympy import symbols, Symbol 
from sympy import Function
from sympy import solve, pretty
from sympy import diff
from sympy import sin, cos, exp, ln, sinh, cosh, sqrt, exp, Abs
from sympy import simplify, trigsimp
from sympy.parsing.sympy_parser import parse_expr


import time

#
#activate if simplify() is not to be used
#
#def simplify(x):
#    return x

def bernoulli_poly(expr):
    """
    Returns a polynom to approximate the Bernoulli function. 
    """
    #return 1.0
    #return 1.0 - expr/2.0
    #return 1.0 - expr/2.0 + expr**2/12.0 - expr**4/720.0 
    return 1.0 - expr/2.0 + expr**2/12.0 - expr**4/720.0 + expr**6/30240.0 
    #return 1.0 - expr/2.0 + expr**2/12.0 - expr**4/720.0 + expr**6/30240.0 - expr**8/1209600.0
    #return expr/(exp(expr) - 1)

def bernoulli_exp(expr):
    """
    Returns the Bernoulli-function.
    """
    return expr/(exp(expr) - 1)


def substituteFunctions(function, sub_functions, left_side, sub_args = None, partial_derivative = None, tabs=0, do_simplify = 0):
    """
    Depending on the value of the argument of a function it is substituted.
    """

    for sub_function in sub_functions:

        found_arguments = []

        for funcs in function.find(Function):

            if(str(type(funcs)) == str(sub_function["function"])):
                    found_arguments.append(funcs.args[0])

        comb = (sub_function, found_arguments, )


        if not found_arguments:
            for j in range(tabs):
                print("\t ", end = '')
            print(left_side, "= ", end = '')

            if  partial_derivative is not None:
                if (do_simplify):
                    print(simplify(diff(function, partial_derivative)).subs(sub_args))
                else:
                    print(diff(function, partial_derivative).subs(sub_args))
            else:
                if (do_simplify):
                    print(simplify(function).subs(sub_args))
                else:
                    print(function.subs(sub_args))
        else:
            for i in range(len(comb[0]["if_states"])**len(comb[1])):
                substitutes = []
                perm = [int(j) for j in str(bin(i))[2:].zfill(len(comb[1]))]
                for j in range(tabs):
                    print("\t ", end = '')
                print("if ", end = '')
                for k in range(len(comb[1])):
                    print('(', end = '')
                    print(comb[0]["if_expr"](comb[1][k]).subs(sub_args), end = '')
                    print(comb[0]["if_states"][perm[k]], end = '')
                    print(comb[0]["value"], end = '')
                    print(')', end = '')
                    if (k < (len(comb[1])-1)):
                        print(" and ", end = '')
                    substitutes.append((comb[0]["function"](comb[1][k]), (comb[0]["subs"][perm[k]](comb[1][k]))), )
                print(":")

                print("###")
                print("### function:", function)
                print("### substitutes:", substitutes)
                print("#### partial_derivative:", partial_derivative)
                print("### sub_args:", sub_args)
                print("###")

                for j in range(tabs+1):
                    print("\t ", end = '')
                print(left_side, "= ", end = '')
                
                if  partial_derivative is not None:
                    part_der = diff(function.subs(substitutes), partial_derivative).subs(sub_args)
                    if (do_simplify):
                        print(simplify(part_der))
                    else:
                        print(part_der)
                else:
                    if (do_simplify):
                        print(simplify(function).subs(substitutes).subs(sub_args))
                    else:
                        print(function.subs(substitutes).subs(sub_args))


def makeUpdate_b(name, functions, search_sub_functions, substitutes):   
    """
    Generates the function vector.
    """
    print("################################")
    print("########### {name} ###########".format(name=name))
    print("################################")

    print("def {name}(Ua, Ub):".format(name=name)) 
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

             substituteFunctions(function[1],
                                 search_sub_functions,
                                 left_side = "parameters.b["+function[0]+"]",
                                 sub_args = substitutes[s],
                                 tabs = 3,
                                 do_simplify = 1
                )

    print() 
    print("\t return None")
    print()
    print()
    print()


def codeGenerator():
    """
    This function generates the python code. 
    """
    bernoulli = Function('bernoulli')

    bernoulli_limit = symbols('parameters.bernoulli_limit')


    Phi_n_k, Phi_n_l = symbols('Phi_n_k, Phi_n_l')
    Phi_p_k, Phi_p_l = symbols('Phi_p_k, Phi_p_l')
    Psi_k, Psi_l     = symbols('Psi_k, Psi_l')

    Chi_k, Chi_l = symbols('Chi_k, Chi_l')
    Eg_k, Eg_l   = symbols('Eg_k, Eg_l')
    Nv_k, Nv_l   = symbols('Nv_k, Nv_l')
    Nc_k, Nc_l   = symbols('Nc_k, Nc_l')

    Psi_m1, Psi_00, Psi_p1       = symbols('parameters.u[3*(i-1)+0], parameters.u[(3*i+0)+0], parameters.u[3*(i+1)+0]')
    Phi_p_m1, Phi_p_00, Phi_p_p1 = symbols('parameters.u[3*(i-1)+1], parameters.u[(3*i+0)+1], parameters.u[3*(i+1)+1]')
    Phi_n_m1, Phi_n_00, Phi_n_p1 = symbols('parameters.u[3*(i-1)+2], parameters.u[(3*i+0)+2], parameters.u[3*(i+1)+2]')


    Chi_m1, Chi_00, Chi_p1    = symbols('parameters.Chi[i-1], parameters.Chi[i+0], parameters.Chi[i+1]')
    Eg_m1, Eg_00, Eg_p1       = symbols('parameters.Eg[i-1], parameters.Eg[i+0], parameters.Eg[i+1]')

    Nv_m1, Nv_00, Nv_p1       = symbols('parameters.Nv[i-1], parameters.Nv[i+0], parameters.Nv[i+1]')
    Nc_m1, Nc_00, Nc_p1       = symbols('parameters.Nc[i-1], parameters.Nc[i+0], parameters.Nc[i+1]')



    dx = symbols('parameters.dx')

    recombination = symbols('parameters.recombination[i]')
    generation    = symbols('parameters.generation[i]')


    C_00 = symbols('parameters.C[i]')
    Epsilon_00 = symbols('parameters.Epsilon[i]')

    np_exp  = Function('np.exp')
    np_sqrt = Function('np.sqrt')
    np_abs  = Function('np.abs')

    kB = symbols('parameters.kB')
    T  = symbols('parameters.T')
    q  = symbols('parameters.q')


    Ut = (kB * T) / q

    mu_n, mu_p = symbols('parameters.mu_n[i], parameters.mu_p[i]')

    Cau = symbols('parameters.Cau[i]')

    ni = Function('ni')
    p = Function('p')
    n = Function('n')



    # ni^2

    ni2 =(Nv_00*Nc_00)*exp(-q*(Eg_00)/(kB*T))
    p  = Nv_00*exp(q*( (-Chi_00 - Eg_00) + Phi_p_00 - Psi_00)/(kB*T))
    n  = Nc_00*exp(q*(        +Chi_00    - Phi_n_00 + Psi_00)/(kB*T))


    #
    # Van Roosbroeck equations
    #

    #poisson = ((Psi_p1 - 2.*Psi_00 + Psi_m1) + q / Epsilon * (C_00 + p -n)*dx**2).subs([(Psi_k, Psi_00),(Psi_l, Psi_p1),(Phi_p_k, Phi_p_00),(Phi_p_l,Phi_p_p1)])
    poisson = ((Psi_p1 - 2.*Psi_00 + Psi_m1) + q / Epsilon_00 * (C_00 + p - n)*dx**2)
    

    # 
    j_p =+(q*mu_p*Ut)/(dx)*(
        +bernoulli(+(Psi_l-Psi_k)/(Ut)) * Nv_k*exp(q*(Phi_p_k - Psi_k + ((-Chi_k - Eg_k) + (-Chi_l - Eg_l))/2.)/(kB*T))
        -bernoulli(-(Psi_l-Psi_k)/(Ut)) * Nv_l*exp(q*(Phi_p_l - Psi_l + ((-Chi_k - Eg_k) + (-Chi_l - Eg_l))/2.)/(kB*T))
    )

    j_n = -(q*mu_n*Ut)/(dx)*(
        +bernoulli(-(Psi_l-Psi_k)/(Ut)) * Nc_k*exp(q*(Psi_k - Phi_n_k - (-Chi_k - Chi_l)/2.)/(kB*T))
        -bernoulli(+(Psi_l-Psi_k)/(Ut)) * Nc_l*exp(q*(Psi_l - Phi_n_l - (-Chi_k - Chi_l)/2.)/(kB*T))
    )


    div_j_p = (
        + j_p.subs([(Psi_k, Psi_00),(Psi_l, Psi_p1),(Phi_p_k, Phi_p_00),(Phi_p_l,Phi_p_p1),(Eg_k, Eg_00),(Eg_l, Eg_p1),(Chi_k, Chi_00),(Chi_l, Chi_p1),(Nv_k, Nv_00),(Nv_l, Nv_p1)])
        - j_p.subs([(Psi_k, Psi_m1),(Psi_l, Psi_00),(Phi_p_k, Phi_p_m1),(Phi_p_l,Phi_p_00),(Eg_k, Eg_m1),(Eg_l, Eg_00),(Chi_k, Chi_m1),(Chi_l, Chi_00),(Nv_k, Nv_m1),(Nv_l, Nv_00)])
        + q*(Cau*(n*p-ni2)-generation)*dx
    )

    div_j_n = (
        + j_n.subs([(Psi_k, Psi_00),(Psi_l, Psi_p1),(Phi_n_k, Phi_n_00),(Phi_n_l,Phi_n_p1),(Chi_k, Chi_00),(Chi_l, Chi_p1),(Nc_k, Nc_00),(Nc_l, Nc_p1)])
        - j_n.subs([(Psi_k, Psi_m1),(Psi_l, Psi_00),(Phi_n_k, Phi_n_m1),(Phi_n_l,Phi_n_00),(Chi_k, Chi_m1),(Chi_l, Chi_00),(Nc_k, Nc_m1),(Nc_l, Nc_00)])
        - q*(Cau*(n*p-ni2)-generation)*dx
    )



    #############################
    # Start of module generation
    #############################

    functions = [
        ["3*i+0",poisson, "poisson"],
        ["3*i+1",div_j_p, "div_j_p"],
        ["3*i+2",div_j_n, "div_j_n"]
    ]


    partial_derivatives = [
        ["3*(i-1)+0",Psi_m1  ],
        ["3*(i+0)+0",Psi_00  ],
        ["3*(i+1)+0",Psi_p1  ],
        ["3*(i-1)+1",Phi_p_m1],
        ["3*(i+0)+1",Phi_p_00],
        ["3*(i+1)+1",Phi_p_p1],
        ["3*(i-1)+2",Phi_n_m1],
        ["3*(i+0)+2",Phi_n_00],
        ["3*(i+1)+2",Phi_n_p1]
    ]



    substitutes = {"left"  : ((Psi_m1, '( ohm_potential(parameters.C[0], parameters.Chi[0], parameters.Eg[0], parameters.Nc[0], parameters.Nv[0]) + Ua)'),
                              (Phi_p_m1, 'Ua'),
                              (Phi_n_m1, 'Ua'),
                              (exp,np_exp),
                              (sqrt,np_sqrt),
                              (Abs,np_abs)
    ),
                   "right" : ((Psi_p1, '( ohm_potential(parameters.C[parameters.n-1], parameters.Chi[parameters.n-1], parameters.Eg[parameters.n-1], parameters.Nc[parameters.n-1], parameters.Nv[parameters.n-1]) + Ub)'
                   ),
                              (Phi_p_p1, 'Ub'),
                              (Phi_n_p1, 'Ub'),
                              (Chi_p1, Chi_00),
                              (Eg_p1, Eg_00),
                              (Nv_p1, Nv_00),
                              (Nc_p1, Nc_00),
                              (exp,np_exp),
                              (sqrt,np_sqrt),
                              (Abs,np_abs)
                   ),
                   "center": ((1,1),
                              (exp,np_exp),
                              (sqrt,np_sqrt),
                              (Abs,np_abs)
                   )
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
    print("import synumses.one_dimension.parameters as parameters")
    print("from   synumses.one_dimension.functions import ohm_potential")

    print()
    print()

    print("#################################")
    print("########### Bernoulli ###########")
    print("#################################")
    print("def bernoulli(x):")

    x = Symbol('x')
    print("\t if (x < ", bernoulli_limit,"):")
    print("\t \t return",bernoulli_poly(x))
    print("\t else:")
    print("\t \t return",bernoulli_exp(x).subs(exp,np_exp))
    print()
    print()

    print("############################################")
    print("########### hole_current_density ###########")
    print("############################################")
    print("def hole_current_density():")
    print()
    print("\t j_p = np.zeros(parameters.n)")
    print()
    print("\t for i in range(0,parameters.n-1):")

    search_sub_function = ({"function" : bernoulli,
                           "if_expr"   : Abs,
                           "if_states" : [' <= ',
                                          ' >  '],
                           "value"     : bernoulli_limit,
                           "subs"      : [bernoulli_poly,
                                          bernoulli_exp]
    },)


    substituteFunctions(j_p,
                        search_sub_function,
                        "j_p[i]",
                        sub_args = 
                        (
                            (Psi_k,   Psi_00), 
                            (Psi_l,   Psi_p1),
                            (Phi_p_k, Phi_p_00),
                            (Phi_p_l, Phi_p_p1),
                            (Chi_k,   Chi_00),
                            (Chi_l,   Chi_p1),
                            (Eg_k,    Eg_00),
                            (Eg_l,    Eg_p1),
                            (Nv_k,    Nv_00),
                            (Nv_l,    Nv_p1),
                            (exp,     np_exp),
                            (Abs,     np_abs)
                        ),
                        tabs = 2,
                        do_simplify = 1
    )

    print()
    print("\t i = parameters.n-1")
    print("\t j_p[i] =  j_p[i-1]")
    print()
    print("\t return j_p")
    print()
    print()

    print("################################################")
    print("########### electron_current_density ###########")
    print("################################################")
    print("def electron_current_density():")
    print()
    print("\t j_n = np.zeros(parameters.n)")
    print()
    print("\t for i in range(0,parameters.n-1):")

    substituteFunctions(j_n,
                        search_sub_function,
                        "j_n[i]",
                        sub_args = 
                        (
                            (Psi_k,   Psi_00), 
                            (Psi_l,   Psi_p1),
                            (Phi_n_k, Phi_n_00),
                            (Phi_n_l, Phi_n_p1),
                            (Chi_k,   Chi_00   ),
                            (Chi_l,   Chi_p1),
                            (Nc_k,    Nc_00),
                            (Nc_l,    Nc_p1),
                            (exp,     np_exp),
                            (Abs,     np_abs)
                        ),
                        tabs = 2,
                        do_simplify = 1
    )



    print()
    print("\t i = parameters.n-1")
    print("\t j_n[i] =  j_n[i-1]")
    print()
    print("\t return j_n")
    print()
    print()


    bernoulli = Function('bernoulli')
    bernoulli_limit = symbols('parameters.bernoulli_limit')

    search_sub_function = {"function"  : bernoulli,
                           "if_expr"   : Abs,
                           "if_states" : [' <= ',
                                          ' >  '],
                           "value"     : bernoulli_limit,
                           "subs"      : [bernoulli_poly,
                                          bernoulli_exp]}

    search_sub_functions = ()
    search_sub_functions = search_sub_functions + (search_sub_function, )

    makeUpdate_b("update_b", functions, search_sub_functions, substitutes)

    
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

        #print("# ",substitutes[s])

        for function in functions:
            print()
            print("\t \t \t #######################")
            print("\t \t \t ### ",function[2]," ###")
            print("\t \t \t #######################")

            for partial_derivative in partial_derivatives:
                if ((s == "left"  and (partial_derivative[1] in [Psi_m1, Phi_p_m1, Phi_n_m1])) or
                    (s == "right" and (partial_derivative[1] in [Psi_p1, Phi_p_p1, Phi_n_p1]))):
                    print("\t \t \t ###"+"parameters.A["
                          +function[0]
                          +","
                          +partial_derivative[0]
                          +"]",
                          " = ... ")
                else:
                    substituteFunctions(-function[1],
                                        search_sub_functions, 
                                        "parameters.A["+function[0]+","+partial_derivative[0]+"]",
                                        substitutes[s],
                                        partial_derivative = partial_derivative[1],
                                        tabs=4,
                                        do_simplify = 1)
                                        

                    



    print()
    print("\t return ")
    print()
    print()

    poisson =  ((Psi_p1 - 2*Psi_00 + Psi_m1) + q / Epsilon_00 * (C_00 + p -n)*dx**2)
    #div_j_p = Phi_p_00
    #div_j_n = Phi_n_00
    div_j_p = Phi_p_00 - Phi_p_m1
    div_j_n = Phi_n_00 - Phi_n_m1

    functions = [
        ["3*i+0",poisson, "poisson"],
        ["3*i+1",div_j_p, "div_j_p"],
        ["3*i+2",div_j_n, "div_j_n"]]

    partial_derivatives = [
        ["3*(i-1)+0",Psi_m1  ],
        ["3*(i+0)+0",Psi_00  ],
        ["3*(i+1)+0",Psi_p1  ],
        ["3*(i-1)+1",Phi_p_m1],
        ["3*(i+0)+1",Phi_p_00],
        ["3*(i+1)+1",Phi_p_p1],
        ["3*(i-1)+2",Phi_n_m1],
        ["3*(i+0)+2",Phi_n_00],
        ["3*(i+1)+2",Phi_n_p1]
    ]


    substitutes = {"left"  : ((Psi_m1, '(ohm_potential(parameters.C[0], parameters.Chi[0], parameters.Eg[0], parameters.Nc[0], parameters.Nv[0]) + Ua)'),
                              (Phi_p_m1, 'Ua'),
                              (Phi_n_m1, 'Ua'),
                              (exp,np_exp),
                              (sqrt,np_sqrt)),
                   "right" : ((Psi_p1, '(ohm_potential(parameters.C[parameters.n-1], parameters.Chi[parameters.n-1], parameters.Eg[parameters.n-1], parameters.Nc[parameters.n-1], parameters.Nv[parameters.n-1]) + Ub)'),
                              (Phi_p_p1, 'Ub'),
                              (Phi_n_p1, 'Ub'),
                              (Eg_p1, Eg_00),
                              (Chi_p1, Chi_00),
                              (Nv_p1, Nv_00),
                              (Nc_p1, Nc_00),
                              (exp,np_exp),
                              (sqrt,np_sqrt)),
                   "center": ((1,1),
                              (exp,np_exp),
                              (sqrt,np_sqrt))
        }



    makeUpdate_b("first_update_b", functions, search_sub_functions, substitutes)


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

        #print("# ",substitutes[s])

        for function in functions:
            print
            print("\t \t \t #######################")
            print("\t \t \t ### ",function[2]," ###")
            print("\t \t \t #######################")
            for partial_derivative in partial_derivatives:
                if ((s == "left"  and (partial_derivative[1] in [Psi_m1, Phi_p_m1, Phi_n_m1])) or
                    (s == "right" and (partial_derivative[1] in [Psi_p1, Phi_p_p1, Phi_n_p1]))):
                    print("\t \t \t ###"+"parameters.A["
                          +function[0]
                          +","
                          +partial_derivative[0]
                          +"]",
                          " = ..."
                    )
                else:
                    substituteFunctions(-function[1],
                                        search_sub_functions, 
                                        "parameters.A["+function[0]+","+partial_derivative[0]+"]",
                                        substitutes[s],
                                        partial_derivative = partial_derivative[1],
                                        tabs=4,
                                        do_simplify = 1)

    print()
    print("\t return ")
    print()
    print()

if __name__ == '__main__':
    codeGenerator()
