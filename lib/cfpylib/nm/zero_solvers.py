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
"""

TEST_MODULE = 0

def secant(f, x0, x1, tol=1.0e-11, limits=[], max_iterations=1000):
    """
    The iterative secant method for zero-finding in one-dimension.
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

def bisection(f, by, uy, tol=1.0e-6):
    """
    The iterative bisection method for zero-finding in one-dimension.
    """
    while abs(uy-by) > tol:
        midpoint = 0.5*(by+uy)
        if f(by) * f(midpoint) > 0:
            by = midpoint
        else:
            uy = midpoint
    return 0.5*(by+uy)

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

        
