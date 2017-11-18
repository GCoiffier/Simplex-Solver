# Optimization & Approximation 2017 : Homework
## Guillaume Coiffier
``` test ```

## How to use the program
I implemented the simplex in python. To run the program, the command is:
`python3 simplex.py [-v] [-d] [-r rule] file `

where options are the following:
- `file` The input file. All inputs file are in the input folder.
- `-v` enables verbose mode : gives a detailed feedback on the execution
- `-d` enables debug mode : gives an even more detailed feedback on the execution
- `-r` rule : choice of rule (default rule is Random). Rules have to be : Random, Bland, MaxCoeff or Custom



To generate random Linear Program, you can use the LPgenerator.py script. To run this script, run the following command:
` python3 LPgenerator.py [-n N] [-m M] [-twophase] [-hollow] outputfile

where options are the following:
- `-n` The number of variables
- `-m` The number of constraints
- `-twophase` Allow the program to generate negative coefficient. This will often result in a 2 phase resolution
- `-hollow` Each coefficient as a 0.5 chance of being zero.

##
