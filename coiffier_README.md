# Optimization & Approximation 2017 : Homework
## Guillaume Coiffier

## How to use the program
I implemented the simplex in python. To run the program, the command is:

    python3 simplex.py [-v] [-d] [-r rule] file

where options are the following:
- `file` The input file. All inputs file are in the input folder.
- `-v` enables verbose mode : gives a detailed feedback on the execution
- `-d` enables debug mode : gives an even more detailed feedback on the execution
- `-r` rule : choice of rule (default rule is Random). Rules have to be : Random, Bland, MaxCoeff or Custom


To generate random Linear Program, you can use the LPgenerator.py script. To run this script, run the following command:

    python3 LPgenerator.py [-n N] [-m M] [-random | -klee-minty D] [-twophase] [-hollow] outputfile

where options are the following:
- `-n` The number of variables
- `-m` The number of constraints
- `-random | -klee-minty D` Generate either a random LP, or the Klee-Minty cube of dimension D.
- `-twophase` In the case of a random generation, allow the program to generate negative coefficient. This will often result in a 2 phase resolution
- `-hollow` In the case of a random generation, each coefficient as a 0.5 chance of being zero.

note that most of the time, the random problems generated with two phases wil be unfeasible or unbounded. An example of a feasible problem is given as coiffier_test_random3.in

## Code architecture
My implementation is divided into 5 python files :
- `simplex.py` The main file, that contain the main() function and the implementation of the simplex method. Main logic of the program is in there.
- `linearProgram.py` The definition of a class representing a linear program
- `tableau.py` The definition of the class Tableau.
- `utilities.py` Utility function to display fractions into the console
- `LPgenerator.py` A small script to help me writing big instances of linear programs in order to generate test files.

## Testing files
test files include:
- Examples from the subject, and some unitary tests
- Some random instances (small and big)
- LP for busDriver, fraction and lumberjack problems
- Klee-Minty's Cube LP to show how some rule can badly perform on some examples
