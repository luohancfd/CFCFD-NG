#!/usr/bin/env python
## \file: bezier_spline.py
## \ingroup: geom

"""
Bezier polynomial functions like those in bezier.c.

P.J. May,June 2004, with RJ Gollan splines, etc

Reference:
    Fujio Yamaguchi
    Curves and Surfaces in Computer Aided Geometric Design
    Springer-Verlag 1988
"""

import sys
from copy import copy
from geom import Vector
from math import fabs,sqrt
from zero_solvers import secant
from bezier import *

class mySpline(object):
    def __init__(self, Bx, By=[], Bz=[]):
        self.Bx = copy(Bx)
        self.By = copy(By)
        self.Bz = copy(Bz)

class myData(object):
    def __init__(self, spline, y_value):
        self.spline = copy(spline)
        self.y_value = y_value

init_guesses = [ [0.25, 0.75], [0.1, 0.3], [0.7, 0.9] ]


def bezier_3_spline(m, p):
    """
    Given m+1 interpolation points, determine the m-segment,
    Bezier polyline that interpolates these points.

    Unlike the C-function, this returns a list of the
    beizer control points used to represent the spline.

    Further comments are in the accompanying C-file.
    Coded in Python: 06-Dec-2004
    """
    
    assert len(p) >= (m+1)

    if m <= 0:
        print 'Cannot attempt to find a Bezier representation'
        print 'with ', m, ' segments.'
        return 'FAIL'

    tolerance = 1.0e-10
        
    # Setup the initial guess for the weight points by using the
    # supplied points. This also sets the end points correctly.

    d = []
    for node in p:
        # Use a Vector object which will not be stored in the Node.nodeList.
        d.append(Vector(node.x, node.y, node.z))
    
    # Apply the Gauss-Seidel interation until the internal
    # weight points converge.

    for j in range(50):
        max_diff = 0.0
        for i in range(1,m):
            old_p = d[i]

            d[i].x = 0.25 * (6.0 * p[i].x - d[i-1].x - d[i+1].x)
            d[i].y = 0.25 * (6.0 * p[i].y - d[i-1].y - d[i+1].y)
            d[i].z = 0.25 * (6.0 * p[i].z - d[i-1].z - d[i+1].z)

            dx = fabs(d[i].x - old_p.x)
            dy = fabs(d[i].y - old_p.y)
            dz = fabs(d[i].z - old_p.z)

            diff = dx + dy + dz
            if ( diff > max_diff ):
                max_diff = diff
        # end for i
        if (max_diff < tolerance):
            break

    # end for j

    if j >= 50:
        print 'Iterations did not converge while using the Guass-Seidel technique'
        print 'to find a bezier representation of a spline.'
        return 'FAIL'

    # Finally, calculate the Bezier segments an pack them away.
    bcp = range(4*m)

    for i in range(m):
        bcp[i*4 + 0] = copy(p[i])

        bcp[i*4 + 1] = Vector()
        bcp[i*4 + 1].x = (2.0 * d[i].x + d[i+1].x) / 3.0
        bcp[i*4 + 1].y = (2.0 * d[i].y + d[i+1].y) / 3.0
        bcp[i*4 + 1].z = (2.0 * d[i].z + d[i+1].z) / 3.0

        bcp[i*4 + 2] = Vector()
        bcp[i*4 + 2].x = (d[i].x + 2.0 * d[i+1].x) / 3.0
        bcp[i*4 + 2].y = (d[i].y + 2.0 * d[i+1].y) / 3.0
        bcp[i*4 + 2].z = (d[i].z + 2.0 * d[i+1].z) / 3.0

        bcp[i*4 + 3] = copy(p[i+1])

    return bcp
    
def spline_eval_based_on_y(spNodes, y_value):
    """
    Given a spline represented by a bezier and a
    y value of interest.  Return the (first) point on
    the spline matches that y value.
    """

    # Find bezier segment of interest
    seg = -1
    n_seg = len(spNodes) / 4
    for i in range(n_seg):
        lower = spNodes[i*4 + 0].y
        upper = spNodes[i*4 + 3].y
        if y_value >= lower and y_value <= upper:
            seg = i
            break

    if seg == -1:
        print 'Did not find y value in the spline interval'
        print 'y= ', y_value
        if y_value < spNodes[0].y:
            # We are below the minimum y
            return ('FAIL', spNodes[0].y, 'FAIL')
        else:
            # We are abov the maximum y_value
            return ('FAIL', spNodes[-1].y, 'FAIL')
        

    b0 = Vector( spNodes[seg*4 + 0].x, spNodes[seg*4 + 0].y, spNodes[seg*4 + 0].z )
    b1 = Vector( spNodes[seg*4 + 1].x, spNodes[seg*4 + 1].y, spNodes[seg*4 + 1].z )
    b2 = Vector( spNodes[seg*4 + 2].x, spNodes[seg*4 + 2].y, spNodes[seg*4 + 2].z )
    b3 = Vector( spNodes[seg*4 + 3].x, spNodes[seg*4 + 3].y, spNodes[seg*4 + 3].z )

    Bx = [ b0.x, b1.x, b2.x, b3.x ]
    By = [ b0.y, b1.y, b2.y, b3.y ]
    Bz = [ b0.z, b1.z, b2.z, b3.z ]

    B = mySpline( Bx, By, Bz )
    funData = myData( B, y_value )


    for guess in init_guesses:
        x0 = guess[0]
        x1 = guess[1]
        t = secant(zero_fun_y, x0, x1, funData, limits=[0.0, 1.0])
        if t != 'FAIL':
            break
            
        
    if t == 'FAIL':
        print 'There is a problem finding the t-parameter'
        sys.exit(-1)

    return bezier3D_eval(t, Bx, By, Bz)

def zero_fun_y(t, funData ):
    """
    The function to find the zero of.
    """

    y_set = funData.y_value
    y_guess = bezier_eval(t, funData.spline.By)

    return (y_set - y_guess)

def bezier2D_length(Bx, By):
    """Return the length of a Bezier curve in the xy-plane."""
    SAMPLES = 1000
    length = 0.0
    dt = 1.0 / SAMPLES
    t = 0.0
    (x_old, y_old, z_old) = bezier3D_eval(t, Bx, By)
    for i in range(1,SAMPLES):
        t = t + dt
        (x_new, y_new, z_new) = bezier3D_eval(t, Bx, By)
        dr = sqrt( (x_new-x_old)*(x_new-x_old) + (y_new-y_old)*(y_new-y_old) )
        length += dr
        x_old = x_new; y_old = y_new; z_old = z_new

    return length

def spline2D_length(bcp):
    """Return the length of a spline in the xy-plane."""
    length = 0.0
    for i in range(0, len(bcp), 4):
        Bx = [ bcp[i].x, bcp[i+1].x, bcp[i+2].x, bcp[i+3].x ]
        By = [ bcp[i].y, bcp[i+1].y, bcp[i+2].y, bcp[i+3].y ]
        seg = bezier2D_length(Bx, By)
        length += seg
    return length

def bezier2D_find_point(length, Bx, By):
    """Find a point on a Bezier curve a given length along."""
    SAMPLES = 1000
    test = 0.0
    dt = 1.0 / SAMPLES
    t = 0.0
    (x_old, y_old, z_old) = bezier3D_eval(t, Bx, By)
    for i in range(1,SAMPLES):
        t = t + dt
        (x_new, y_new, z_new) = bezier3D_eval(t, Bx, By)
        dr = sqrt( (x_new-x_old)*(x_new-x_old) + (y_new-y_old)*(y_new-y_old) )
        test += dr
        if test >= length:
            return Node(x_new, y_new)
        x_old = x_new; y_old = y_new;
        
    return 'FAIL'
    
def spline2D_find_point(length, bcp):
    """Find a point on a spline a given length along."""
    test = 0.0
    for i in range(0, len(bcp), 4):
        Bx = [ bcp[i].x, bcp[i+1].x, bcp[i+2].x, bcp[i+3].x ]
        By = [ bcp[i].y, bcp[i+1].y, bcp[i+2].y, bcp[i+3].y ]
        seg = bezier2D_length(Bx, By)
        test += seg
        if test >= length:
            # We have found the segment
            bezLength = seg - (test - length)
            return bezier2D_find_point(bezLength, Bx, By)

    return 'FAIL'

# -------------------------------------------------------------------

if __name__ == '__main__':
    print "Bezier demo begin..."
    print bezier3D_eval(0.66667, [0.0, 2.0, 4.0, 6.0],
                        [-2.0, 1.0, 4.0, 7.0])

    print "Testing Bezier length eval..."
    Bx = [0.0, 2.0, 4.0, 6.0]
    By = [0.0, 8.0/3.0, 16.0/3.0, 8.0]
    print bezier2D_length(Bx, By)
    

    print "Done."
    
