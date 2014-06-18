# adapti.py
"""
Adaptive quadrature using Newton-Cotes 5- and 3-point rules.

.. Author: PA Jacobs, School of Engineering, UQ
.. MECH2700 demonstration code.
.. Version 18-Nov-2003
.. Version 26-jan-2013 Sphinx docs
.. Version 14-Jun-2014 halve tolerances when passing them down

Example transcript::

    peterj@laval ~/cfcfd3/lib/cfpylib/nm $ python adapti.py 
    Begin adapti...
    Estimates of pi/4:  0.785398154922 0.785398173551
    errors: 8.47578540686e-09 -1.01539638919e-08
    number of function calls: 315 55
    Done.
"""

def rinteg(f, a, b, tol):
    """
    Apply Newton-Cotes 5- and 3-point quadrature rules to the segment [a,b].

    :param f: user-supplied function, f(x)
    :param a, b: range of integration
    :param tol: maximum difference between rules,
                above which the range is split
    :returns: integral of f(x) from a to b.
    """
    dx = b - a
    f0 = f(a)
    f1 = f(a + 0.25 * dx)
    f2 = f(a + 0.5 * dx)
    f3 = f(a + 0.75 * dx)
    f4 = f(b)
    I2 = dx/6.0 * (f0 + 4 * f2 + f4)
    I4 = dx/90.0 * (7*f0 + 32*f1 + 12*f2 + 32*f3 + 7*f4)
    if abs(I4 - I2) > tol:
        mid = 0.5 * (a + b)
        I = rinteg(f, a, mid, tol/2.0) + rinteg(f, mid, b, tol/2.0)
    else:
        I = I4
    return I

from math import pi, sqrt

if __name__ == "__main__":
    count1 = 0
    def fun1(x):
        global count1
        count1 += 1
        if abs(x) < 1.0:
            return sqrt(1.0 - x * x)
        else:
            return 0.0

    count2 = 0
    def fun2(x):
        global count2
        count2 += 1
        return 1.0 / (1.0 + x * x)

    print "Begin adapti..."
    a = 0.0; b = 1.0;
    pi4_1 = rinteg(fun1, a, b, 1.0e-5)
    pi4_2 = rinteg(fun2, a, b, 1.0e-5)
    print "Estimates of pi/4: ", pi4_1, pi4_2
    print "errors:", pi/4 - pi4_1, pi/4 - pi4_2
    print "number of function calls:", count1, count2

    print "Done."
    
