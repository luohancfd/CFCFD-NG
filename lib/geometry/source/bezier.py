#!/usr/bin/env python
## \file: bezier.py
## \ingroup: geom

"""
Bezier polynomial functions like those in bezier.c.

P.J. May,June 2004.

Reference:
    Fujio Yamaguchi
    Curves and Surfaces in Computer Aided Geometric Design
    Springer-Verlag 1988
"""

from copy import copy

def bezier_eval(t, B):
    """
    Evaluate the Bezier polynomial at parametric loation 0<=t<=1.

    For an n-th order curve, we expect n+1 points.
    The de Casteljau algorithm is used for simplicity.
    """
    assert len(B) > 1
    assert t > -0.1 and t < 1.01

    n = len(B) - 1   # order of the polynomial
    Q = copy(B)      # work array

    # Now, generate one new level at a time;
    # over-writing the work array at each level.
    for k in range(n):
        for i in range(n-k):
            Q[i] = (1.0 - t) * Q[i] + t * Q[i+1]

    return Q[0]

def bezier3D_eval(t, Bx, By=[], Bz=[]):
    """
    Evaluate the Bezier polynomial for 1, 2 or 3D space.
    """
    assert len(Bx) > 1
    n = len(Bx) - 1  # order of the polynomial
    if not By: By = [0.0] * (n+1)
    if not Bz: Bz = [0.0] * (n+1)
    assert len(By) == n+1
    assert len(Bz) == n+1
    return (bezier_eval(t,Bx), bezier_eval(t,By), bezier_eval(t,Bz))

def bezier_add_one_point(B):
    """
    Adds one control point and returns the new list of control points.

    The algorithm in Section 5.1.8 on page 208 of Yamaguchi is used.
    This effectively increases the order of the polynomial by one.
    """
    n = len(B) - 1
    assert n >= 1
    Q = [B[0], ]  # First point is a special case
    for i in range(1,n+1):
        Qnew = (i * B[i-1] + (n + 1 - i) * B[i]) / (n + 1.0)
        Q.append(Qnew)
    Q.append(B[n]) # Last point is also a special case
    return Q

# -------------------------------------------------------------------

if __name__ == '__main__':
    print "Bezier demo begin..."
    print bezier3D_eval(0.66667, [0.0, 2.0, 4.0, 6.0],
                        [-2.0, 1.0, 4.0, 7.0])
    print "Done."
    
