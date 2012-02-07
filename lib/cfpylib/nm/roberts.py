#! /bin/env python
"""
roberts.py: Node distribution and coordinate stretching functions.

These functions should behave the same as the C code functions.

Author PA Jacobs

Version: 1.0, 22-Mar-2005
"""

try:
    from numpy import arange, power
except:
    try:
        from Numeric import arange, power
    except:
        print "Could import neither numpy nor Numeric."

def roberts(eta, alpha, beta):
    """
    Computes the stretched coordinate in the range [0.0..1.0]
    using the boundary-layer-like transformation devised by Roberts.
    
    :param eta: unstretched coordinate, 0 <= eta <= 1.0
    :param beta: stretching factor (more stretching as beta --> 1.0)
    :param alpha: location of stretching:

        | alpha = 0.5: clustering of nodes at both extremes of eta
        | alpha = 0.0: nodes will be clustered near eta=1.0

    Works for both scalars and arrays.
    """
    lmbda = (beta + 1.0) / (beta - 1.0)
    # Note that power is a ufunc from Numeric
    lmbda = power(lmbda, ((eta - alpha)/(1.0 - alpha)))
    etabar = (beta + 2.0 * alpha) * lmbda - beta + 2.0 * alpha
    etabar = etabar / ((2.0 * alpha + 1.0) * (1.0 + lmbda))
    return etabar


def distribute_points_1(t1, t2, n, end1, end2, beta):
    """
    Generate a set of points nonuniformly distributed from t1 to t2.

    :param t1: parameter value  1
    :param t2: parameter value  2
    :param n: number of intervals (the number of points is n+1)
    :param end1: (int) cluster flag for end 1:

        | ==0 points are not clustered to end 1
        | ==1 points ARE clustered to end 1

    :param end2: (int) cluster flag for end 2:

        | ==0 points are not clustered to end 2
        | ==1 points ARE clustered to end 2

    :param beta: grid stretching parameter:

        | 1 < beta < +inf : points are clustered 
        | The closer to 1, the more the clustering.
        | beta < 1 for no clustering.
 
    :returns: t[0:n] an array of distributed values.
    """ 
    # Decide on stretching parameters for Robert's transform.
    alpha = 0.5;
    reverse = 0;
    cluster = 1;
    if (end1 == 0 and end2 == 0) or beta < 1.0:
        cluster = 0
    if end1 == 1 and end2 == 1:
        alpha = 0.5
    if end1 == 1 and end2 == 0: 
        reverse = 1
        alpha = 0.0
    if end1 == 0 and end2 == 1:
        reverse = 0
        alpha = 0.0
    # Compute the grid points as an array.
    # The intermediate parameter is uniformly distributed.
    del_eta = 1.0 / n
    i = arange(n+1)
    eta = del_eta * i;
    # Cluster the points.
    if cluster:
        if reverse: eta = 1.0 - eta;
        etabar = roberts(eta, alpha, beta);
        if reverse: etabar = 1.0 - etabar;
    else:
        etabar = eta;
    # Compute the parameter value within the given end-points.
    x = (1.0 - etabar) * t1 + etabar * t2;
    return x

#------------------------------------------------------------------

if __name__ == "__main__":
    print "Begin roberts demo..."
    a = 0.5
    b = 1.1
    print "Scalar use..."
    for eta in arange(0.0, 1.0, 0.1):
        print "eta=", eta, "roberts=", roberts(eta, a, b)
    print "Vector use..."
    eta = arange(0.0, 1.0, 0.1)
    print "eta=", eta, "roberts=", roberts(eta, a, b)
    print "Distribute points..."
    print "x=", distribute_points_1(0.0, 1.0, 5, 1, 0, 1.1)
    print "x=", distribute_points_1(0.0, 1.0, 5, 0, 1, 1.1)
    print "x=", distribute_points_1(0.0, 1.0, 5, 1, 1, 1.1)
    print "x=", distribute_points_1(0.0, 1.0, 5, 0, 0, 1.1)
    print "Done."
