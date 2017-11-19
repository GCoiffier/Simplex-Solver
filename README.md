# Optimization & Approximation 2017 : Homework
## Guillaume Coiffier

## How to use the program
I implemented the simplex in python. To run the program, the command is:
    `python3 simplex.py [-v] [-d] [-r rule] file `

where options are the following:
- `file` The input file. All inputs file are in the input folder.
- `-v` enables verbose mode : gives a detailed feedback on the execution
- `-d` enables debug mode : gives an even more detailed feedback on the execution
- `-r` rule : choice of rule (default rule is Random). Rules have to be : Random, Bland, MaxCoeff or Custom


To generate random Linear Program, you can use the LPgenerator.py script. To run this script, run the following command:

    `python3 LPgenerator.py [-n N] [-m M] [-random | -klee-minty D] [-twophase] [-hollow] outputfile`

where options are the following:
- `-n` The number of variables
- `-m` The number of constraints
- `-random | -klee-minty D` Generate either a random LP, or the Klee-Minty cube of dimension D.
- `-twophase` In the case of a random generation, allow the program to generate negative coefficient. This will often result in a 2 phase resolution
- `-hollow` In the case of a random generation, each coefficient as a 0.5 chance of being zero.

## Custom pivot rule
The custom pivot rule I choosed is the "steepest edge" rule.


## General Remarks
#### Minimization problem
This simplex implementation only works for maximisation problems.
To get a minimisation problem, you will have to maximize the opposite of the objective function you had.
Be careful, because by doing this, you will get at the end the opposite of the optimal value !

#### Artifical variables problem
When doing the phaseI-phaseII method, there exists cases where, at the end of phase I, there remains some
artificial variables in the basis. If those variables are not equal to zero, then we know that the LP
is unfeasible. However, they can be equal to zero. In that case, we need to perform some additionnal pivot
to get the variables out of the basis (and if no pivot can be performed, one can show that we can suppress the
corresponding constraint without changing the solution of the LP).

This method has been implemented in my simplex, but I did not find any example where artificial variables are
present in the basis at the end of phase I, so I was unable to test it.
