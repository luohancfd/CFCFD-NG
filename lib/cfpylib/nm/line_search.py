# line_search.py
"""
Implementation of an algorithm for optimization from Gerald and Wheatley.

A class demo for mech2700.
PJ, 23-Oct-2008 
"""

def minimize(f, a, b, tolerance=1.0e-4):
    """
    Returns the bracket xL,xR containing the minimum of the function f.

    :param f: a user supplied function
    :param a,b: the original bracket containing a minimum
    :param tolerance: the final size of the bracket.
                      It should not be set too small.
    """
    r = 0.618034
    xL = a + (1-r)*(b-a)
    xR = a + r*(b-a)
    FL = f(xL)
    FR = f(xR)

    while (xR - xL) > tolerance:
        if FR > FL:
            b = xR
            xR = xL
            FR = FL
            xL = a + (1-r)*(b-a)
            FL = f(xL)
        else:
            a = xL
            xL = xR
            FL = FR
            xR = a + r*(b-a)
            FR = f(xR)
        # print "xL=%g  xR=%g  FL=%g  FR=%g" % (xL, xR, FL, FR)
    return xL, xR


if __name__ == '__main__':
    print "Begin..."
    from math import exp, cos

    def fdemo(x):
        return exp(x) + 2 - cos(x)

    x = -1.4721
    print "x=", x, "f(", x, ")=", fdemo(x)

    a = -3
    b = 1
    print "bracket=", minimize(fdemo, a, b, 1.0e-6)

    print "Done."
