"""
ridder.py
Solve a nonlinear equation f(x)=0 using the method of Ridder.

Peter J. 
mech3750 demo code 12-Mar-2014
"""

from math import sqrt, copysign

def solve(f, x1, x2, tol=1.0e-9):
    """
    Locate a root of f(x) by subdividing the original range,
    assuming a linear model of the underlying transformed function.

    Input:
        f: user-supplied function f(x)
        x1: first end of range
        x2: other end of range
    Returns:
        x, a point near the root.
    """
    assert callable(f), "User-supplied function must be callable."
    f1 = f(x1); f2 = f(x2)
    if abs(f1) == 0.0: return x1
    if abs(f2) == 0.0: return x2
    assert x1 != x2, "Bad initial range given to bracket."
    assert f1 * f2 < 0.0, "Range does not clearly bracket a root."
    while abs(x2 - x1) > tol:
        x3 = 0.5*(x1+x2)
        f3 = f(x3)
        if f3 == 0.0: return x3
        eps = (f3 + copysign(sqrt(f3*f3-f1*f2),f2))/f2
        x4 = x3 - f3*eps*(x1-x3)/(f1 - eps*f3)
        f4 = f(x4)
        if f4 == 0.0: return x4
        # Contract the bracket.
        if f3*f2 < 0.0:
            if f4*f2 < 0.0:
                x1 = x4; f1 = f4
            else:
                x1 = x3; f1 = f3
                x2 = x4; f2 = f4
        else:
            if f4*f1 < 0.0:
                x2 = x4; f2 = f4
            else:
                x1 = x4; f1 = f4
                x2 = x3; f2 = f3
    return x4

if __name__ == '__main__':
    print "Begin self-test of Ridder's method..."
    #
    from math import pow, sin, exp
    def test_fun_1(x):
        return ( pow(x,3) + pow(x,2) - 3*x - 3 )
    print ''
    print 'Example from Gerald and Wheatley, p. 45'
    print 'Solve f(x) = x^3 + x^2 - 3x -3 = 0 with initial'
    print 'guesses of x0 = 1 and x1 = 2.'
    print 'Final result x = ', solve(test_fun_1, 1, 2)
    print 'Gerald and Wheatley report x = 1.732051'
    print ''
    #
    def test_fun_2(x):
        return ( 3*x + sin(x) - exp(x) )
    print 'Example from Gerald and Wheatley, p.45 also'
    print 'Solve f(x) = 3*x + sin(x) - e^x = 0 with initial'
    print 'guesses of x0 = 0 and x1 = 1.'
    print 'Final result x = ', solve(test_fun_2, 0, 1)
    print 'Gerald and Wheatley report x = 0.3604217'
    print "Done."
