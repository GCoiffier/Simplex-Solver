from fractions import *
import numpy as np
from .utilities import *

# ======================= LinearProgram Class ==================================
class LinearProgram:
    """
    Inner representation of a linear program.
    Contains the following datas :
        - nbVar : the initial number of variables
        - nbConst : number of constraints
        - objectiveFunction : the vector c such that transpose(c)*x is the objective value
        - objectiveVector : the vector b such that constraint i is <= b
        - constraintMatrix
    """

    # _____ Parsing and initialisation _____
    def __init__(self,filename="test1.in"):
        # retrieving data from file
        f = open(filename, 'r')
        input_data = [ s.strip() for s in f.readlines()]
        f.close()

        # formatting to build a LinearProgram object
        self.nbVar = int(input_data[0])
        self.nbConst = int(input_data[1])
        self.objectiveFunction = np.array([convert(u) for u in input_data[2].strip().split()])
        self.constraintVector = np.array([convert(u) for u in input_data[3].strip().split()])
        self.constraintMatrix = np.array([Fraction(0,1)]*self.nbVar*self.nbConst).reshape((self.nbConst,self.nbVar))
        for i in range(self.nbConst):
            l = [convert(u) for u in input_data[4+i].strip().split()]
            for j in range(self.nbVar):
                self.constraintMatrix[i,j]=l[j]

        self.need_2_phases = False # Real value set at computation of the standard form of the LP
        f.close()

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
                    first_coef = False
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
                    first_coef = False
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
