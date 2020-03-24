
from sympy import symbols, Symbol 
from sympy import Function
from sympy import solve, pretty
from sympy import diff
from sympy import sin, cos, exp, ln, sinh, cosh, sqrt, exp, Abs
from sympy import simplify, trigsimp
from sympy.parsing.sympy_parser import parse_expr
import itertools



def substituteFunctions(function, sub_function, left_side, sub_args = [], partial_derivative = None, tabs=0, do_simplify = 0):
    """
    Depending on the value of the argument of a function it is substituted.
    """


    # Search for all arguments of (sub_function["function"]

    comb = []

    
    found_arguments = []
    for funcs in function.find(Function): #loop all functions

        if(str(type(funcs)) == str(sub_function["function"])):
                found_arguments.append(funcs.args[0])

    comb = (sub_function, found_arguments)

    print("len(found_arguments):", len(found_arguments))
    print("len(sub_function()", len(sub_function["if_states"])+1)

    print(list(range(len(sub_function["if_states"]))))

    perm = list(itertools.product(range(len(sub_function["if_states"])+1), repeat = len(found_arguments)))

    

    if not found_arguments:

        for ii in range(tabs):
            print("\t", end = '')

        print(left_side, "= ", end = '')
        print(function)

    else:
        
        for i in range(len(perm)):

            for ii in range(tabs):
                print("\t", end = '')

            if (i==0):
                print("if (", end = '')
            elif(i<(len(perm)-1)):
                print("elif (", end = '')
            else:
                print("else", end = '')

            stateSet = 0;
            substitutes = []
            for j in range(len(found_arguments)): #number of arguments

                #print("comb[0][]:", comb[0]["function"](comb[1][j]))
                #print("comb[0][]:", comb[0]["if_subs"][perm[i][j]](comb[1][j]))

                #substitutes.append((comb[0]["function"](comb[1][j]), comb[0]["if_subs"][perm[i][j]](comb[1][j])))

                if (j==0):
                    if (perm[i][j] != len(sub_function["if_states"])):
                        #print("(perm:",perm[i][j], "arg:", j, ")",end = '')
                        print("(", comb[1][j], comb[0]["if_states"][perm[i][j]], comb[0]["if_values"][perm[i][j]], ")", end = '')
                        substitutes.append((comb[0]["function"](comb[1][j]), comb[0]["if_subs"][perm[i][j]](comb[1][j])))
                        stateSet=1
                    else:
                        substitutes.append((comb[0]["function"](comb[1][j]), comb[0]["else_sub"](comb[1][j])))
                        
                else:
                    if (perm[i][j] != len(sub_function["if_states"])):
                        substitutes.append((comb[0]["function"](comb[1][j]), comb[0]["if_subs"][perm[i][j]](comb[1][j])))
                        if (stateSet==1):
                            #print(" AND ", "(perm:", perm[i][j], "arg:", j, ")",end = '')
                            print(" AND (", comb[1][j], comb[0]["if_states"][perm[i][j]], comb[0]["if_values"][perm[i][j]], ")", end = '')
                            
                        else:
                           #print("(perm:", perm[i][j], "arg:", j, ")",end = '') 
                           print("(", comb[1][j], comb[0]["if_states"][perm[i][j]], comb[0]["if_values"][perm[i][j]], ")", end = '')
                           stateSet = 1;
                    else:
                       substitutes.append((comb[0]["function"](comb[1][j]), comb[0]["else_sub"](comb[1][j]))) 
                
            if(i<(len(perm)-1)):
                print("):\n", end = '')
            else:
                print(":\n", end = '')

            for ii in range(tabs):
                print("\t", end = '')

            print("\t" , end = '')
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


            print()



if __name__ == '__main__':

    x,y,z = symbols('x, y, z')
    k,l,m = symbols('k, l, m')

    bernoulli_limit = symbols('parameters.bernoulli_limit')

    f11 = Function('f11')
    f12 = Function('f12')
    f13 = Function('f13')

    f1 = Function('f1')
    f2 = Function('f2')
    f3 = Function('f3')
    f4 = Function('f4')
    f5 = Function('f5')


   
    f_ges = f11(x) + f12(y) + f11(z) + f11(y)


    search_sub_function = ({"function"  : f11,
                            "if_states" : ['<', '<'],
                            "if_values" : [-bernoulli_limit, bernoulli_limit],
                            "if_subs"   : [f1, f2],
                            "else_sub"  : f1
    })


    substituteFunctions(f_ges, search_sub_function, "A", tabs=1)

    #codeGenerator()
