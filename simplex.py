################################################################################
#                      Simplex algorithm Implementation                        #
#                    Guillaume Coiffier, M1IF @ENS de Lyon                     #
#                      Optimisation & Approximation Course                     #
#                                  2017                                        #
################################################################################

from fractions import *
import numpy as np
import argparse
from random import randint

pivotRules = {"Random", "Bland", "MaxCoeff", "Custom"} # Implemented pivot rules

# ========== Exception Definitions =============================================

class Unbounded(Exception):
    pass

class Infeasible(Exception):
    pass

# ========== Utility function over fractions ===================================
def convert(u):
    """
    u is a string representing a rationnal number.
    Outputs the fraction object whose value is described in u
    """
    u = u.split("/")
    if len(u)==1 : #integer
        return Fraction(int(u[0]),1)
    else :
        return Fraction(int(u[0]), int(u[1]))

def frac_print(u):
    if u.denominator == 1:
        return str(u.numerator)
    else :
        return str(u.numerator) +"/"+ str(u.denominator)


# ======================= LinearProgram Class ==================================
class LinearProgram:

    # _____ Initialisation and input_____
    def __init__(self,n,m,c,b,A):
        self.need_2_phases = False
        self.nbVar = n
        self.nbConst = m
        self.objectiveFunction = c
        self.constraintVector = b
        self.constraintMatrix = A

    # _____ Output fontions ______
    def __str__(self):
        output_string = "OUTPUT \nThe input linear program is: \n\n"

        obj_func_line = "Maximize  "
        first_coef = True
        for i in range(len(self.objectiveFunction)):
            if (self.objectiveFunction[i]!=0):
                if (not first_coef and self.objectiveFunction[i]>0):
                    obj_func_line += "+ {0}x_{1} ".format(self.objectiveFunction[i], i+1)
                elif (first_coef and self.objectiveFunction[i]>0):
                    first_coef = False
                    obj_func_line += "{0}x_{1} ".format(self.objectiveFunction[i], i+1)
                elif (self.objectiveFunction[i]<0):
                    obj_func_line += "{0}x_{1} ".format(self.objectiveFunction[i], i+1)
        output_string+= obj_func_line+"\n"

        const_line = "Such that "
        for j in range(len(self.constraintMatrix)):
            first_coef = True
            for i in range(len(self.constraintMatrix[j,:])):
                if (not first_coef and self.constraintMatrix[j,i]>0):
                    const_line += "+ {0}x_{1} ".format(frac_print(self.constraintMatrix[j,i]), i+1)
                elif (first_coef and self.constraintMatrix[j,i]>0):
                    first_coef = False
                    const_line += "{0}x_{1} ".format(frac_print(self.constraintMatrix[j,i]), i+1)
                elif (self.constraintMatrix[j,i]<0):
                    const_line += "{0}x_{1} ".format(frac_print(self.constraintMatrix[j,i]), i+1)

            const_line += " <= " + frac_print(self.constraintVector[j])
            const_line += "\n          "
        output_string+= const_line
        non_neg_line = ", ".join(["x_{}".format(i+1) for i in range(self.nbVar)])
        non_neg_line += " are non-negative\n"
        output_string += non_neg_line
        return output_string

    def __repr__(self):
        return self.__str__()


# ============================== Tableau Class =================================
class Tableau:

    def __init__(self, lp):
        """
        Builds the initial tableau of the linear program 'lp'
        This implies tranforming the LP from canonic to standard from.
        This constructor also add the artificial variables in order to run phase 1 of the simplex
        """
        self.nbPivot = 0

        top_row = [Fraction(0,1)]*(lp.nbVar+lp.nbConst)
        n = lp.nbVar+lp.nbConst+1 # will be the total number of columns in the tableau

        self.nonBasicVariables = []
        self.basicVariables = []
        self.varAssocToConstraint = [-1]*(lp.nbConst+1) # gives the basic variable expressed by constraint i+1

        artificialConstRows = []

        for ind,x in enumerate(lp.constraintVector):
            if x<0:
                # we add an artificial variable to run phase 1
                lp.need_2_phases = True # We will need two phases to run the simplex
                top_row.append(-1)
                self.basicVariables.append(n) # the artificial variable created is basic
                self.varAssocToConstraint[ind+1]=n
                artificialConstRows.append(ind+1)
                n +=1
            else:
                slack_var = ind+lp.nbVar+1
                self.varAssocToConstraint[ind+1]=slack_var
                self.basicVariables.append(slack_var) # the slack variable is basic

        self.nonBasicVariables = [x for x in range(1,n) if x not in self.basicVariables]
        top_row.append(0) # initial value of the objective function

        self.width = n # width of tableau = number of columns
        self.height = lp.nbConst+1 # height of tableau = number of rows

        self.tab = np.array([Fraction(0,1)]*self.width*self.height).reshape((self.height, self.width))
        self.tab[0] = top_row

        artificalVarCount = 0
        for i in range(1,self.height): #constraint matrix
            for j in range(lp.nbVar):
                self.tab[i,j] = lp.constraintMatrix[i-1,j]
            self.tab[i,j+i] = Fraction(1,1)
            b = lp.constraintVector[i-1]
            self.tab[i,-1] = b
            if b<0: # one artificial variable is associated with this constraint
                self.tab[i] *= -1 # we want only >0 numbers in the right hand side
                self.tab[i,lp.nbConst+lp.nbVar+artificalVarCount] = 1
                artificalVarCount +=1

        if not lp.need_2_phases: #top line
            self.tab[0] = np.concatenate([lp.objectiveFunction, np.array( [Fraction(0,1)]*(self.width-lp.nbVar))])
        else:
            for line in artificialConstRows:
                self.tab[0,:] += self.tab[line,:]

    def reload_utility_function(self,objfunc):
        """
        Changes the utility function of the Tableau.
        Useful in the transition from phase 1 to phase 2
        """
        self.tab[0] = np.concatenate([objfunc, np.array( [Fraction(0,1)]*(self.width-len(objfunc)))])
        if debugMode:
            print("Import initial objective function :")
            print(self)
        for x in self.basicVariables:
            if self.tab[0,x-1]!=0:
                i = np.argmax(self.tab[1::,x-1])+1 #the line with the 1
                self.tab[0] -= self.tab[i]*self.tab[0,x-1]/self.tab[i,x-1]

    def get_basic(self):
        return self.basicVariables

    def print_basic(self):
        """ For outputing the final base """
        return ", ".join(self.basicVariables)

    def get_non_basic(self):
        return self.nonBasicVariables

    def get_value_of_solution(self):
        """ The current value of the function to maximize """
        return -self.tab[0,self.width-1]

    def get_solution_variables(self, n):
        """ Returns a string containing the value of the n first variables and their values in the current tableau"""
        l = []
        for i in range(n):
            if (i+1) in self.basicVariables:
                l.append("x_{0} = {1}".format(i+1, self.tab[self.varAssocToConstraint[i+1],-1]))
            else:
                l.append("x_{0} = 0".format(i+1))
        return ", ".join(l)

    def do_pivot(self, enteringVar, leavingVar, leavingInd):
        """
        Apply the pivot.
        enteringVar -> the variable that will replace leavingVar in the basis.
        leavingInd -> the index of the constraint representing leavingVar in the tableau
        """
        self.nbPivot += 1
        self.varAssocToConstraint[leavingInd]=enteringVar

        # update the basicVariables and nonBasicVariables lists
        self.basicVariables.remove(leavingVar)
        self.nonBasicVariables.append(leavingVar)
        self.basicVariables.append(enteringVar)
        self.nonBasicVariables.remove(enteringVar)

        # do pivot on the matrix
        self.tab[leavingInd,:] /= self.tab[leavingInd,enteringVar-1] # renormalize
        for i in range(self.height): #apply pivot
            if i != leavingInd:
                self.tab[i,:] -= self.tab[leavingInd,:]*self.tab[i,enteringVar-1]


    def __str__(self):
        # determine spacing between columns
        spacing=0
        for j in range(self.width):
            for i in range(self.height):
                x = frac_print(self.tab[i,j])
                spacing = max(spacing, len(x)+2)

        output_string = ""

        def print_row(i):
            row = ""
            for j in range(self.width-1):
                x = frac_print(self.tab[i,j])
                row += x+" "*(spacing-len(x))
            row+= " |  "
            row += frac_print(self.tab[i,self.width-1])+"\n"
            return row

        top_row = print_row(0) # top row
        output_string += top_row
        output_string += ("-"*(len(top_row)+2))+"\n" # separator
        for i in range(1,self.height):# other rows
            output_string += print_row(i)
        return output_string

    def __repr__(self):
        return self.__str__()

# ========== Simplex functions =================================================

def simplex_choose_entering(tb):
    """ depends on the pivot rule """
    n = -1
    if rule=="Random":
        non_neg = [x for x in tb.get_non_basic() if tb.tab[0,(x-1)]>0]
        if non_neg :
            n = non_neg[randint(0,len(non_neg)-1)]
    elif rule=="Bland":
        for x in tb.get_non_basic():
            if tb.tab[0,(x-1)]>0:
                n=x
                break;
    elif rule=="MaxCoeff":
        n = np.argmax(tb.tab[0,0:-1])+1
        if tb.tab[0,n-1]<=0:
            n=-1
    elif rule=="Custom":
        raise NotImplementedError
    else:
        raise Exception("Pivot rule is not valid !")
    if verboseMode and n!=-1:
        print("The entering variable is x_{0}".format(n))
    return n

def simplex_choose_leaving(tab, enteringVar):
    n = -1
    m = float("inf")
    for j in range(1,tab.height):
        if tab.tab[j,enteringVar-1]>0:
            maxval = tab.tab[j,-1]/tab.tab[j,enteringVar-1]
            if 0<=maxval<m:
                m=maxval
                n=j
    if n==-1: # no upper bound
        raise Unbounded
    var = tab.varAssocToConstraint[n] # find the basic variable associated with row n
    if verboseMode:
        print("The leaving variable is x_{0} \n".format(var))
    return var,n

def simplex_one_phase(tab):
    while True:
        if debugMode:
            print("Basic variables : " + str(tab.get_basic()))
            print("Non basic variables : "+ str(tab.get_non_basic()))
            print("Variables associated to constraints : " + str(tab.varAssocToConstraint[1::]) +"\n")
        inVar = simplex_choose_entering(tab)
        if inVar==-1:
            # We are done : no variable can improve the solution
            return tab
        else:
            outVar,outInd = simplex_choose_leaving(tab, inVar)
            tab.do_pivot(inVar, outVar, outInd)
        if verboseMode:
            print(tab)
    return


def simplex_solve(lp):
    print(lp)
    tab = Tableau(lp)
    if verboseMode:
        print("The initial tableau is : \n")
        print(tab)
    try:
        if lp.need_2_phases:
            # compute phase 1
            if verboseMode:
                print("=========== PHASE 1 ==========\n")
            tab = simplex_one_phase(tab)
            val = tab.get_value_of_solution()
            if (val!=0):
                raise Infeasible
            tab.reload_utility_function(lp.objectiveFunction)
            if verboseMode:
                print("\n========== PHASE 2 ==========\n")
                print(tab)
        else:
            print("The point (0,...,0) is a feasible solution. Only one phase is needed")
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

# ========== Parsing function ============================

def lp_parse(filename):
    """ opens and read a .in file that represent a linear program """

    # retrieving data from file
    f = open(filename, 'r')
    input_data = [ s.strip() for s in f.readlines()]
    f.close()

    # formatting to build a LinearProgram object
    n = int(input_data[0])
    m = int(input_data[1])
    c = np.array([convert(u) for u in input_data[2].strip().split()])
    b = np.array([convert(u) for u in input_data[3].strip().split()])
    A = np.array([Fraction(0,1)]*n*m).reshape((m,n))
    for i in range(m):
        l = [convert(u) for u in input_data[4+i].strip().split()]
        for j in range(n):
            A[i,j]=l[j]
    return LinearProgram(n,m,c,b,A)

# ================== MAIN ======================================================
verboseMode = False
debugMode = False
rule = "Random"

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
    my_lp = lp_parse(filename) # open and parse the lp
    simplex_solve(my_lp) # solve the lp using simplex algorithms
