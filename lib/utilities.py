from fractions import *

# ========== Utility function over fractions ===================================
def convert(u):
    """
    u is a string representing a rationnal number (of general form [-]a[/b])
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
