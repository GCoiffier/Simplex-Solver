################################################################################
#                      Simplex algorithm Implementation                        #
#                    Guillaume Coiffier, M1IF @ENS de Lyon                     #
#                      Optimisation & Approximation Course                     #
#                                  2017                                        #
################################################################################

import argparse
from time import *
from random import randint
from lib import *
import numpy as np

verboseMode = False
debugMode = False
rule = "Random"

pivotRules = {"Random", "Bland", "MaxCoeff", "Custom"}

# ========== Exception Definitions =============================================

class Unbounded(Exception):
    pass

class Infeasible(Exception):
    pass

# =============== Simplex algorithm functions ==================================

def simplex_choose_entering(tab):
    """ depends on the pivot rule """
    n = -1
    if rule=="Random":
        non_neg = [x for x in tab.get_non_basic() if tab[0,(x-1)]>0]
        if non_neg :
            n = non_neg[randint(0,len(non_neg)-1)]
    elif rule=="Bland":
        for x in tab.get_non_basic():
            if tab[0,(x-1)]>0:
                n=x
                break;
    elif rule=="MaxCoeff":
        n = np.argmax(tab[0,0:-1])+1
        if tab[0,n-1]<=0:
            n=-1
    elif rule=="Custom":
        t = tab[0,0:-1].copy()
        for i in range(len(t)):
            if np.dot(tab[:,i], tab[:,i])!=0:
                t[i] /= np.dot(tab[:,i], tab[:,i])
            else:
                t[i]=-1
        n = np.argmax(t)+1
        if tab[0,n-1]<=0:
            n=-1
    else:
        raise Exception("Pivot rule is not valid !")
    if verboseMode and n!=-1:
        print("The entering variable is x_{0}".format(n))
    return n

def simplex_choose_leaving(tab, enteringVar):
    n = -1
    m = float("inf")
    for j in range(1,tab.height):
        if tab[j,enteringVar-1]>0:
            maxval = tab[j,-1]/tab[j,enteringVar-1]
            if 0<=maxval<m:
                m=maxval
                n=j
    if n==-1: # no upper bound
        raise Unbounded
    var = tab.varAssocToConstraint[n] # find the basic variable associated with row n
    assert(tab.constraintAssocToVar[var]==n) # just to make sure
    if verboseMode:
        print("The leaving variable is x_{0} \n".format(var))
    return var

def simplex_one_phase(tab):
    while True:
        if debugMode:
            print("Basic variables : " + str(tab.get_basic()))
            print("Non basic variables : "+ str(tab.get_non_basic()))
            print("Variables associated to constraints : " + str(tab.varAssocToConstraint) +"\n")
        inVar = simplex_choose_entering(tab)
        if inVar==-1:
            # We are done : no variable can improve the solution
            return tab
        else:
            outVar = simplex_choose_leaving(tab, inVar)
            tab.do_pivot(inVar, outVar)
        if verboseMode:
            print(tab)
    return


def simplex_solve(lp):
    start_time = time()
    print(lp)
    tab = Tableau(lp)
    if verboseMode:
        print("The initial tableau is : \n")
        print(tab)
        if debugMode:
            print("Artificial variables : "+str(tab.artificialVariables)+"\n")
    try:
        if lp.need_2_phases:
            # compute phase 1
            if verboseMode:
                print("=========== PHASE 1 ==========\n")
            tab = simplex_one_phase(tab)
            val = tab.get_value_of_solution()
            if (val!=0):
                raise Infeasible
            if debugMode:
                print("END OF PHASE 1\n")
            tab.transition_phaseI_phaseII(lp.objectiveFunction, verboseMode, debugMode)
            if verboseMode:
                print("\n========== PHASE 2 ==========\n")
                print("After reloading the initial objective function, the tableau is:\n")
                print(tab)
        else:
            print("The point (0,...,0) is a feasible solution. Only one phase is needed\n")
        # compute phase 2
        tab = simplex_one_phase(tab)

    except Infeasible:
        print("This linear program in INFEASIBLE")
        return
    except Unbounded:
        print("This linear program is UNBOUNDED")
        return

    print("An optimal solution is : {0}".format(tab.get_solution_variables(lp.nbVar)))
    print("The value of the objective for this solution is : {0}".format(tab.get_value_of_solution()) )
    print("The number of pivots is : {0}".format(tab.nbPivot))
    print("The pivot rule used : {0}".format(rule))
    print("The calculation took {0:.3f} seconds".format(time()-start_time))

# ================== MAIN ======================================================
if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description='Implementation of the simplex algorithm. Done by Guillaume Coiffier. M1IF Opt&Approx 2017-2018 @ENS de Lyon')
    argparser.add_argument('filename', help="name of the source file.")
    argparser.add_argument('-rule', help="specify the pivot's rule. Default is random", default="Random")
    argparser.add_argument('-v', action="store_true", help="enables verbose mode")
    argparser.add_argument('-d', action="store_true", help="enables debug mode")

    options=argparser.parse_args()

    if options.v:
        verboseMode = True

    if options.d:
        verboseMode = True
        debugMode = True

    if options.rule:
        if options.rule not in pivotRules:
            print("The rule '{0}' does not refer to any implemented rule. \n Possible rules are {1} \n".format(options.rule, ", ".join(pivotRules)))
            raise Exception("No correct rule specified. Program will stop")
        else:
            rule = options.rule
            if debugMode:
                print("The following pivot rule will be used : {}".format(rule))

    filename = options.filename
    my_lp = LinearProgram(filename) # open and parse the lp
    simplex_solve(my_lp) # solve the lp using simplex algorithms
