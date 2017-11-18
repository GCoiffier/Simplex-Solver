""" Generator of random LP for testing """

import argparse
from random import randint

def generate_lp(output_file, nbVar, nbConst, twophase, hollow):
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

# ================== MAIN ======================================================
if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description='Linear Program generator. Done by Guillaume Coiffier. M1IF Opt&Approx 2017-2018 @ENS de Lyon')
    argparser.add_argument('filename', help="name of the output thefile.")
    argparser.add_argument('-n', help="number of variables")
    argparser.add_argument('-m', help="number of constraints")
    argparser.add_argument('-twophase', action="store_true", help="Allow the generator to output a LP that need 2 phases")
    argparser.add_argument('-hollow', action="store_true", help="Create a matrix with a lot of 0s")
    options=argparser.parse_args()

    filename = options.filename

    generate_lp(filename, int(options.n), int(options.m), options.twophase, options.hollow)
