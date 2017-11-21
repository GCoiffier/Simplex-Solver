## Custom pivot rule
The custom pivot rule I choosed is the "steepest edge" rule. From empirical tests, this rule seems to run a bit faster than Bland's rule. However, like the maxCoeff rule, we have no garantee on the termination (the rule might cycle infinitely).
Here is a sum up of the number of pivots of each rule on some examples :

                                Bland       MaxCoeff (Dantzig)      Custom (Steepest Edge)        Random

coiffier_test_random3.in          42                49                      50                    random
coiffier_test_random2.in          115               20                      29                    random
coiffier_test_random.in           2                 2                       1                     random
coiffier_test_random4.in          8                 7                       9                     random
coiffier_partiel.in               3                 3                       3                     random
coiffier_klee_minty_10.in         177               1023                    1                     random
coiffier_klee_minty_20.in         21891             1048575                 1                     random
coiffier_busDriver.in             8                 8                       13                    random


## Minimization problem
This simplex implementation only works for maximisation problems.
To get a minimisation problem, you will have to maximize the opposite of the objective function you had.
Be careful, because by doing this, you will get at the end the opposite of the optimal value !

## The Klee-Minty Cube
The Klee-Minty Cube LP is a problem designed to show that the Dantzig's pivot rule (max coefficient) can perform an exponential number of pivots.
For the cube of dimension D (that has D variables and D constraints), we obtain the following numbers of pivots :

#### D=10
- Dantzig (Max Coeff) : 1023 (2^D -1)
- Bland : 177
- Steepest Edge (Custom) : 1

#### D=20 (Warning : Very long computation time : ~1h for Dantzig)
- Dantzig (Max Coeff) : 1048575
- Bland : 21891
- Steepest Edge (Custom) : 1

It is very interesting to see that the Steepest Edge pivot rule finds the optimal point with only one pivot.


## Artifical variables problem
When doing the phaseI-phaseII method, there exists cases where, at the end of phase I, there remains some
artificial variables in the basis. If those variables are not equal to zero, then we know that the LP
is unfeasible. However, they can be equal to zero. In that case, we need to perform some additionnal pivot
to get the variables out of the basis (pivot between the artificial variable and the associated slack variable)

This is the case of the "fraction" problem.
NOTE: There might be cases where the slack variable is already in the basis. In those cases, it's all out of my hands...
