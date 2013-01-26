#!/usr/bin/env python
"""
zero_solvers.py

.. Author: Rowan J Gollan

.. Versions:
   06-Dec-2004
   08-May-2011: Dan's bisection_method added by PJ
   16-Apr-2012: PJ, make more efficient by not evaluating f redundantly
                Also, make the code more compact (so that the full 
                function fits in the editor window).

Example transcript::

    python ~/e3bin/cfpylib/nm/zero_solvers.py
    Begin zero_solvers self-test...

    Test function 1.
    ----------------
    Example from Gerald and Wheatley, p. 45
    Solve f(x) = x^3 + x^2 - 3x -3 = 0 with initial
    guesses of x0 = 1 and x1 = 2.
    Begin function call secant()...

    Iteration 	 x0 		x1 		x2 	 f(x2) 
    -----------------------------------------------------------------------
    1 	  1.000000 	 2.000000 	 1.571429 	 -1.364431e+00
    2 	  2.000000 	 1.571429 	 1.705411 	 -2.477451e-01
    3 	  1.571429 	 1.705411 	 1.735136 	 2.925540e-02
    4 	  1.705411 	 1.735136 	 1.731996 	 -5.151769e-04
    5 	  1.735136 	 1.731996 	 1.732051 	 -1.039000e-06
    6 	  1.731996 	 1.732051 	 1.732051 	 3.702993e-11
    7 	  1.732051 	 1.732051 	 1.732051 	 1.776357e-15
    -----------------------------------------------------------------------
    Final result x =  1.73205080757
    Gerald and Wheatley report x = 1.732051
    Using bisection... x = 1.73205080757

    Test function 2.
    ----------------
    Example from Gerald and Wheatley, p.45
    Solve f(x) = 3*x + sin(x) - e^x = 0 with initial
    guesses of x0 = 0 and x1 = 1.
    Begin function call secant()...

    Iteration 	 x0 		x1 		x2 	 f(x2) 
    -----------------------------------------------------------------------
    1 	  1.000000 	 0.000000 	 0.470990 	 2.651588e-01
    2 	  0.000000 	 0.470990 	 0.372277 	 2.953367e-02
    3 	  0.470990 	 0.372277 	 0.359904 	 -1.294813e-03
    4 	  0.372277 	 0.359904 	 0.360424 	 5.530053e-06
    5 	  0.359904 	 0.360424 	 0.360422 	 1.021329e-09
    6 	  0.360424 	 0.360422 	 0.360422 	 -8.881784e-16
    -----------------------------------------------------------------------
    Final result x =  0.36042170296
    Gerald and Wheatley report x = 0.3604217
    Using bisection... x = 0.360421702961
"""

TEST_MODULE = 0

def secant(f, x0, x1, tol=1.0e-11, limits=[], max_iterations=1000):
    """
    The iterative secant method for zero-finding in one-dimension.

    :param f: user-defined function f(x)
    :param x0: first guess
    :param x1: second guess, presumably close to x0
    :tol: stopping tolerance for f(x)=0
    :max_iterations: to stop the iterations running forever, just in case...

    :returns: x such that f(x)=0 or the string 'FAIL'
    """
    # We're going to arrange x0 as the oldest (furtherest) point
    # and x1 and the closer-to-the-solution point.
    # x2, when we compute it, will be the newest sample point.
    f0 = f(x0); f1 = f(x1)
    if abs(f0) < abs(f1):
        x0, f0, x1, f1 = x1, f1, x0, f0
    for i in range(max_iterations):
        try :
            x2 = x1 - f1 * (x0 - x1) / (f0 - f1)
        except ZeroDivisionError:
            return 'FAIL'
        if limits != []:
            x2 = max(limits[0], x2)
            x2 = min(limits[1], x2)
        f2 = f(x2)
        if TEST_MODULE == 1:
            print '  %d \t  %f \t %f \t %f \t %e' % (i+1, x0, x1, x2, f2 )
        x0, f0, x1, f1 = x1, f1, x2, f2
        if abs(f2) < tol: return x2
    if TEST_MODULE == 1:
        print 'The secant method did not converge after ', i+1, ' iterations'
    return 'FAIL'

# -------------------------------------------------------------------

def bisection(f, bx, ux, tol=1.0e-6):
    """
    The iterative bisection method for zero-finding in one-dimension.

    :param f: user-defined function f(x)
    :param bx: bottom-limit of bracket
    :param ux: upper-limit of bracket
    :tol: stopping tolerance on bracket size

    :returns: x such that f(x)=0
    """
    while abs(ux-bx) > tol:
        midpoint = 0.5*(bx+ux)
        if f(bx) * f(midpoint) > 0:
            bx = midpoint
        else:
            ux = midpoint
    return 0.5*(bx+ux)

# -------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin zero_solvers self-test..."

    from math import pow, sin, exp
    def test_fun_1(x):
        return ( pow(x,3) + pow(x,2) - 3*x - 3 )
    print ''
    print 'Test function 1.'
    print '----------------'
    print 'Example from Gerald and Wheatley, p. 45'
    print 'Solve f(x) = x^3 + x^2 - 3x -3 = 0 with initial'
    print 'guesses of x0 = 1 and x1 = 2.'
    print 'Begin function call secant()...'
    print ''
    print 'Iteration \t x0 \t\tx1 \t\tx2 \t f(x2) '
    print '-----------------------------------------------------------------------'
    TEST_MODULE = 1
    x2 = secant(test_fun_1, 1, 2)
    print '-----------------------------------------------------------------------'
    print 'Final result x = ',x2
    print 'Gerald and Wheatley report x = 1.732051'
    print 'Using bisection... x =', bisection(test_fun_1, 1.0, 2.0, tol=1.0e-11)
    print ''

    def test_fun_2(x):
        return ( 3*x + sin(x) - exp(x) )
    print 'Test function 2.'
    print '----------------'
    print 'Example from Gerald and Wheatley, p.45'
    print 'Solve f(x) = 3*x + sin(x) - e^x = 0 with initial'
    print 'guesses of x0 = 0 and x1 = 1.'
    print 'Begin function call secant()...'
    print ''
    print 'Iteration \t x0 \t\tx1 \t\tx2 \t f(x2) '
    print '-----------------------------------------------------------------------'
    x2 = secant(test_fun_2, 0, 1)
    print '-----------------------------------------------------------------------'
    print 'Final result x = ',x2
    print 'Gerald and Wheatley report x = 0.3604217'
    print 'Using bisection... x =', bisection(test_fun_2, 0.0, 1.0, tol=1.0e-11)
    print ''

        
