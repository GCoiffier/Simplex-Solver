""" Generator of random LP for testing """

import argparse
from random import randint

def generate_random(output_file, nbVar, nbConst, twophase, hollow):
    thefile = open(output_file, 'w')
    thefile.write("{0}\n".format(nbVar))
    thefile.write("{0}\n".format(nbConst))
    if twophase:
        mini = -100
    else:
        mini = 0
    maxi = 100

    # Generate objective function
    for i in range(nbVar):
        if hollow and randint(0,1)==1:
            thefile.write("0 ")
        else:
            thefile.write(str(randint(mini,maxi))+" ")

    thefile.write("\n")

    # Generate constraint vector
    for i in range(nbConst):
        thefile.write(str(randint(mini,maxi))+" ")
    thefile.write("\n")

    # Generate matrix
    for i in range(nbConst):
        for j in range(nbVar):
            if hollow and randint(0,1)==1:
                thefile.write("0 ")
            else:
                thefile.write(str(randint(mini,maxi))+" ")
        thefile.write("\n")

    thefile.close()

def generate_klee_minty(output_file,d):
        """ see https://en.wikipedia.org/wiki/Klee%E2%80%93Minty_cube """
        thefile = open(output_file, 'w')
        thefile.write(str(d)+"\n")
        thefile.write(str(d)+"\n")

        multipleof5 = [5]*(d+1)
        multipleof2 = [1]*(d+1)
        for i in range(1,d+1):
            multipleof5[i]= multipleof5[i-1]*5
            multipleof2[i]= multipleof2[i-1]*2
        multipleof5 = [ str(x) for x in multipleof5 ]
        multipleof2 = [ str(x) for x in multipleof2 ]
        thefile.write(" ".join([multipleof2[d-1-i] for i in range(d)])) # Objective function
        thefile.write("\n")
        thefile.write(" ".join([multipleof5[i] for i in range(d)])) # Constraint vector
        thefile.write("\n")

        multipleof2[1]="1"
        for i in range(1,d+1): # Constraint matrix
            l = multipleof2[i:0:-1]
            for j in range(d-i):
                l.append("0")
            thefile.write(" ".join(l))
            thefile.write("\n")

        thefile.close()

# ================== MAIN ======================================================
if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description='Linear Program generator. Done by Guillaume Coiffier. M1IF Opt&Approx 2017-2018 @ENS de Lyon')
    argparser.add_argument('filename', help="name of the output thefile.")
    argparser.add_argument('-n', help="number of variables")
    argparser.add_argument('-m', help="number of constraints")
    argparser.add_argument('-random', action="store_true", help="generate a random LP")
    argparser.add_argument('-twophase', action="store_true", help="Random generation parameter. Allow the generator to output a LP that need 2 phases")
    argparser.add_argument('-hollow', action="store_true", help="Random generation parameter. Create a matrix with a lot of 0s")
    argparser.add_argument('-klee-minty', help="Generate the Klee Minty cube of dimension d")
    options=argparser.parse_args()

    filename = options.filename
    if options.random:
        generate_random(filename, int(options.n), int(options.m), options.twophase, options.hollow)
    elif options.klee_minty is not None:
        generate_klee_minty(filename,int(options.klee_minty))
