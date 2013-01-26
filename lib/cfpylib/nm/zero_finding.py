#! /usr/bin/env python
"""
zero_finding.py

.. Author:  Brendan O'Flaherty

.. Versions:
    23 Oct 2007 -- File created
    12 Nov 2007 -- Implemented bisection and golden methods

Example transcript::

    $ python ~/e3bin/cfpylib/nm/zero_finding.py
    bisection method
    find x for f0(x) < 1.000000e-06 for x(2.20:4.08)
    f0(3.141593) = 4.494084e-07

    golden section method
    find x for min(f1(x)) < 1.000000e-06 for x(0.94:2.20)
    f1(1.570796) = -1.000000
"""

import sys

from math import pi, sin, tan
#from pylab import figure, plot, savefig

# ---------------------------------------------------------
# Numerical derivatives.

def derivatives(f, xval, dx, args=()):
    y0 = f(xval, args);
    y1 = f(xval+dx, args);
    y2 = f(xval-dx, args);
    
    dydx = (y1 - y2)/(2*dx);
    #ddydxx = (y1 - 2.0*y0 + y2)/(dx*dx); 

    return y0, dydx#, ddydxx

# ----------------------------------------------------------
# Contains a number of constrained zero finding methods.

def bisection_root(f, x0, x1, tol, args=(), max_iter=100):
    """
    Applies the Bisection-Root method to find f(x) = 0.
    """
    assert x1 > x0

    y0 = f(x0, args)
    y1 = f(x1, args)
    
    if (y0*y1) >= 0.0:
        print "error in bisection root:"
        print "x-bounds [%g, %g] correspond to y-bounds [%g, %g]," % (x0, x1, y0, y1)
        print "which do not fall either side of zero."
        sys.exit()

    niter = 0
    #print "we have found x-bounds [%g, %g] and y-bounds of [%g, %g]" % (x0, x1, y0, y1)
    while abs(y1-y0) > tol:
        xval = 0.5*(x0+x1)
        yval = f(xval, args)
        #print "we have found the midpoint x=%g, y=%g" % (xval, yval)
        if (y0*yval) < 0.0:
            #print "new range [%g, %g]" % (x0, xval)
	    x1 = xval
            y1 = yval
	else:
            #print "new range [%g, %g]" % (xval, x1)
	    x0 = xval
            y0 = yval
        niter = niter + 1
        if (niter >= max_iter):
            print "# maximum number of interations reached!"
            print "# breaking..."
            break

    return 0.5*(x0+x1)

def muller_root(f, x0, x1, tol, args=(), max_iter=100):
    """
    Applies the Bisection-Root method to find f(x) = 0.
    """
    assert x1 > x0

    y0 = f(x0, args)
    y1 = f(x1, args)
    
    if (y0*y1) >= 0.0:
        print "error in bisection root:"
        print "x-bounds [%g, %g] correspond to y-bounds [%g, %g]," % (x0, x1, y0, y1)
        print "which do not fall either side of zero."
        sys.exit()

    niter = 0
    #print "we have found x-bounds [%g, %g] and y-bounds of [%g, %g]" % (x0, x1, y0, y1)
    y_root = 1.0e3
    
    while abs(y_root) > tol:
        xval = 0.5*(x0+x1)
        yval = f(xval, args)

        q = (x1 - xval)/(xval - x0)
        A = q*y1 - q*(1 + q)*yval + q*q*y0
        B = (2*q + 1)*y1 - (1 + q)*(1 + q)*yval + q*q*y0
        C = (1 + q)*y1

        den =  B*B - 4*A*C
        if (den < 0):
            print("error in muller root:")
            print(" complex root when evaluating polynomial.")
            print(" exiting to system...")
            sys.exit()

        x_new1 = x1 - (x1 - xval)*(2*C/(B + den**0.5))
        x_new2 = x1 - (x1 - xval)*(2*C/(B - den**0.5))

        # now we have five points
        # the original bounds, the midpoint, and the two polynomial roots
        # first check which root is within the original bounds

        if ((x_new1 > x0) and (x_new1 < x1)):
            x_root = x_new1
        else:
            x_root = x_new2
        
        # now we have four points
        # the original bounds, the midpoint, and a good approximation.
        # find the new bounds.

        y_root = f(x_root, args)

        if (y0*yval <= 0.0): # use left interval
            x1 = xval
            y1 = yval
            if (y_root*y1 <= 0.0): # use right interval
                x0 = x_root
                y0 = y_root
            else: # use the left interval
                x1 = x_root
                y1 = y_root
            
        else: # use right interval
            x0 = xval
            y0 = yval
            if (y_root*y1 <= 0.0): # use right interval
                x0 = x_root
                y0 = y_root
            else: # use left interval
                x1 = x_root
                y1 = y_root
            
        niter = niter + 1
        if (niter >= max_iter):
            print "# maximum number of interations reached!"
            print "# breaking..."
            break

    return x_root

def newton_root(f, x0, x1, tol, args=()):
    """
    Applies Newton's method to find f(x) = 0.
    """
    assert x1 > x0
    
    y0 = f(x0, args)
    y1 = f(x1, args)

    if (y0*y1) >= 0.0:
        print "error in Newton root"
        print "x-bounds must fall either side of root."
        sys.exit()

    converged = False
    x = (x0+x1)/2.0
    for i in range(50):
        y, dydx = derivatives(f, x, tol, args)
        x -= y/dydx
        if abs(y) < tol:
            converged = True
            break
        
    if converged:
        return x
    else:
        print "Newtons method failed to converge after 50 iterations"
        sys.exit()
        
# ----------------------------------------------------------
# Contains a number of constrained minima finding methods.
# These are effectively zero-finding methods for the derivative of a function f(x)

def golden_section(f, a, c, tol, args=()):
    """
    Applies the Golden-Section search to find min(f(x)).
    """
    g = (3.0 - 5**0.5)/2.0
    b1 = a + g*(c - a)
    b2 = a + (1 - g)*(c - a)
    f1 = f(b1)
    f2 = f(b2)

    while (c - a) > tol*(a + c):
        if f1 <= f2:
            c = b2
            b2 =  b1
            f2 = f1
            b1 = a + g*(c - a)
            f1 = f(b1)
        else:
            a = b1
            b1 = b2
            f1 = f2
            b2 = a + (1.0 - g)*(c - a)
            f2 = f(b2)

    return (a + c)/2.0

# ----------------------------------------------------------

def f0(x, args=()):
    return tan(x)

def test0():
    tol = 1.0e-6
    xmin = 0.7*pi
    xmax = 1.3*pi

    print "bisection method"
    print "find x for f0(x) < %e for x(%3.2f:%3.2f)" % (tol, xmin, xmax)
    x0 = bisection_root(f0, xmin, xmax, tol)
    print "f0(%f) = %e" % (x0, f0(x0))
    print
    
def f1(x, args=()):
    return -sin(x)

def test1():
    tol = 1.0e-6
    xmin = 0.3*pi
    xmax = 0.7*pi

    print "golden section method"
    print "find x for min(f1(x)) < %e for x(%3.2f:%3.2f)" % (tol, xmin, xmax)
    x0 = golden_section(f1, xmin, xmax, tol)
    print "f1(%f) = %f" % (x0, f1(x0))
    print
    
if __name__ == '__main__':
    test0()
    test1()
    
