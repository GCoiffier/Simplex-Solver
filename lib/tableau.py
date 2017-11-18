import numpy as np
from fractions import *
from .utilities import *

# ============================== Tableau Class =================================
class Tableau:
    """
    Implementation of the full tableau method in the simplex algorithm.
    The tableau of a Linear Program is a concise way to represent its data.
    Simplex algorithm then do gaussian pivots on this tableau.
    """

    def __init__(self, lp):
        """
        Builds the initial tableau of the linear program 'lp'
        This implies tranforming the LP from canonic to standard from.
        This constructor also add artificial variables in order to run phase 1
        of the simplex when it is necessary
        """
        self.nbPivot = 0 # Counter for output

        top_row = [Fraction(0,1)]*(lp.nbVar+lp.nbConst)
        n = lp.nbVar+lp.nbConst+1 # will be the total number of columns in the tableau

        self.nonBasicVariables = []
        self.basicVariables = []
        self.artificialVariables = []
        self.varAssocToConstraint = [-1]*(lp.nbConst+1) # gives the basic variable expressed by constraint i+1

        artificialConstRows = []

        for ind,x in enumerate(lp.constraintVector):
            if x<0:
                # we add an artificial variable to run phase 1
                lp.need_2_phases = True # We will need two phases to run the simplex
                top_row.append(-1)
                self.artificialVariables.append(n)
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

        self.data = np.array([Fraction(0,1)]*self.width*self.height).reshape((self.height, self.width))
        self.data[0] = top_row

        artificalVarCount = 0
        for i in range(1,self.height): #constraint matrix
            for j in range(lp.nbVar):
                self.data[i,j] = lp.constraintMatrix[i-1,j]
            self.data[i,j+i] = Fraction(1,1)
            b = lp.constraintVector[i-1]
            self.data[i,-1] = b
            if b<0: # one artificial variable is associated with this constraint
                self.data[i] *= -1 # we want only >0 numbers in the right hand side
                self.data[i,lp.nbConst+lp.nbVar+artificalVarCount] = 1
                artificalVarCount +=1
        # top line
        if not lp.need_2_phases:
            self.data[0] = np.concatenate([lp.objectiveFunction, np.array( [Fraction(0,1)]*(self.width-lp.nbVar))])
        else:
            for line in artificialConstRows:
                self.data[0,:] += self.data[line,:]

    def transition_phaseI_phaseII(self, objfunc, verboseMode, debugMode):
        """
        Changes the utility function of the Tableau.
        Useful in the transition from phase 1 to phase 2
        """
        self.data[0] = np.concatenate([objfunc, np.array( [Fraction(0,1)]*(self.width-len(objfunc)))])
        if debugMode:
            print("Import initial objective function :")
            print(self)
        for x in self.basicVariables:
            if self.data[0,x-1]!=0:
                i = np.argmax(self.data[1::,x-1])+1 #the line with the 1
                self.data[0] -= self.data[i]*self.data[0,x-1]/self.data[i,x-1]

    def get_basic(self):
        return self.basicVariables

    def print_basic(self):
        """ For outputing the final base """
        return ", ".join(self.basicVariables)

    def get_non_basic(self):
        return self.nonBasicVariables

    def get_value_of_solution(self):
        """ The current value of the function to maximize """
        return -self.data[0,self.width-1]

    def get_solution_variables(self, n):
        """ Returns a string containing the value of the n first variables and their values in the current tableau"""
        l = []
        for i in range(n):
            if (i+1) in self.basicVariables:
                l.append("x_{0} = {1}".format(i+1, self.data[self.varAssocToConstraint[i+1],-1]))
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
        self.data[leavingInd,:] /= self.data[leavingInd,enteringVar-1] # renormalize
        for i in range(self.height): #apply pivot
            if i != leavingInd:
                self.data[i,:] -= self.data[leavingInd,:]*self.data[i,enteringVar-1]


    def __str__(self):
        # determine spacing between columns
        spacing=0
        for j in range(self.width):
            for i in range(self.height):
                x = frac_print(self.data[i,j])
                spacing = max(spacing, len(x)+2)

        output_string = ""

        def print_row(i):
            row = ""
            for j in range(self.width-1):
                x = frac_print(self.data[i,j])
                row += x+" "*(spacing-len(x))
            row+= " |  "
            row += frac_print(self.data[i,self.width-1])+"\n"
            return row

        top_row = print_row(0) # top row
        output_string += top_row
        output_string += ("-"*(len(top_row)+2))+"\n" # separator
        for i in range(1,self.height):# other rows
            output_string += print_row(i)
        return output_string

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value
