#!/usr/bin/env python
"""
zero_solvers.py

Author: Rowan J Gollan

Versions:

   | 06-Dec-2004
   | 08-May-2011: Dan's bisection_method added by PJ
"""

from copy import copy
from math import fabs, pow, sin, exp

TEST_MODULE = 0
ZERO_VAL = 1.0e-6
max_iterations = 1000


def secant(f, x0, x1, tol=1.0e-11, limits=[]):
    """
    The iterative secant method for zero-finding in one-dimension.
    """
    f0 = f(x0)
    f1 = f(x1)

    if fabs(f0) < fabs(f1):
        temp = x0; x0 = x1; x1 = temp;

    for i in range(max_iterations):
        f0 = f(x0); f1 = f(x1);
        try :
            x2 = x1 - f1 * (x0 - x1) / (f0 - f1)
        except ZeroDivisionError:
            return 'FAIL'
            
        if limits != []:
            if x2 < limits[0]:
                x2 = limits[0]
            if x2 > limits[1]:
                x2 = limits[1]
        
        if TEST_MODULE == 1:
            print '  %d \t  %f \t %f \t %f \t %e' % (i+1, x0, x1, x2, f(x2) )
        x0 = x1; x1 = x2

        if fabs( f(x2) ) < tol:
            return x2

    if i >= (max_iterations-1):
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

def test_fun_1(x):
    return ( pow(x,3) + pow(x,2) - 3*x - 3 )

def test_fun_2(x):
    return ( 3*x + sin(x) - exp(x) )

# -------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin zero_solvers self-test..."
    print ''
    print 'Test function 1.'
    print '----------------'
    print 'Example from Gerald and Wheatley, p. 45'
    print 'Solve f(x) = x^3 + x^2 - 3x -3 = 0 with initial'
    print 'guesses of x0 = 1 and x1 = 2.'
    print 'Begin function call secant()...'
    print ''
    print 'Iteration \t x0 \t\tx1 \t\tx2 \t F(x2) '
    print '-----------------------------------------------------------------------'
    TEST_MODULE = 1
    x2 = secant(test_fun_1, 1, 2)
    print '-----------------------------------------------------------------------'
    print 'Final result x = ',x2
    print 'Gerald and Wheatley report x = 1.732051'
    print 'Using bisection... x =', bisection(test_fun_1, 1.0, 2.0, tol=1.0e-11)
    print ''

    print 'Test function 2.'
    print '----------------'
    print 'Example from Gerald and Wheatley, p.45'
    print 'Solve f(x) = 3*x + sin(x) - e^x = 0 with initial'
    print 'guesses of x0 = 0 and x1 = 1.'
    print 'Begin function call secant()...'
    print ''
    print 'Iteration \t x0 \t\tx1 \t\tx2 \t F(x2) '
    print '-----------------------------------------------------------------------'
    x2 = secant(test_fun_2, 0, 1)
    print '-----------------------------------------------------------------------'
    print 'Final result x = ',x2
    print 'Gerald and Wheatley report x = 0.3604217'
    print 'Using bisection... x =', bisection(test_fun_2, 0.0, 1.0, tol=1.0e-11)
    print ''

        
