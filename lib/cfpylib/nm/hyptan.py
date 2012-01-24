#!/bin/env python
## \file hyptan.py
## \brief Python script demonstrating hyperbolic tangent grid clustering
## \author Adriaan Window
## \version November 2005
##
"""

This script introduces the idea of hyperbolic tangents as grid clustering
methods of for finite volume flow calculations.
The benefit of hyperbolic tangent clustering is that it allows the user to 
specify the desired spacing as a dimensioned quantity as opposed to an
arbitrary scaling parameter beta. This clustering method also allows differing
values of clustering at each end of the line providing extra flexibility over
the basic symmetry clustering already provided in mb_cns.


The following steps are performed in generating the clustered node 
distribution:
1. User specifies the spacing desired at either end of the line and the number
of nodes to be generated
2. Given the length of the line from the geometry of the block, the spacing
values are normalised against the length.
3. Given that different spacing at each end is desired the line is divided into
two equal lengths and thus equal number of nodes.
4. For each half, a scaling factor is found which is applied to the function
    tanh( beta * x) where x is the vector along which nodes are to be generated
5. Depending on the user spacing compared to uniform spacing, the two halves
are merged into one array of distributed nodes
6. This combined array is then smoothed in the middle 20% of the line to 
eliminate any discontinuities in the distribution due to the initial splitting


1. Ideally I would like this code to also control clustering about some 
arbitrary midpoint. Intuitively, this would involve splitting the nodes about 
that point rather than simply halving. However, depending on the location this
may leave a sparse set of nodes on one interval causing extreme discontinuities
and perhaps even failing at the desired clustering with the given algorithm.
2. Currently the half split occasionally results in spacing larger than 
desired, particularly on the end of the interval where the largest spacing
has be specified. Ultimately this end should have less nodes than the tightly
spaced end, and so the half splitting of nodes is not sufficient. Some better
method for splitting the nodes along the interval is required, based partly on 
the spacing specified by the user (ie a ratio or such).
3. Following the integration of this method into the mb_cns gridding system, 
a set of checkpoints should be put in place to ensure that neighbouring blocks
have the same spaceing in adjacent cells and that the gradients are not 
steep. This could otherwise induce numerical diffusion.

VERSION HISTORY

2005-11-24  -- Working version of clustering completed, still some things to 
sort out as per discussion above
2005-11-29  -- Polished hyptan.py for integration with mb_cns and scriptit.py
2005-12-09  -- Improved ability of algorithm to deal with larger specified
spacing at either end of the interval. This will also allow for internal 
clustering once implemented.
2006-01-11  -- Fixed issue where intervals would cluster incorrectly due to
bad combinations of interval length and number of nodes. Now, clustering is
forced to be similar at both ends (for the time being).
"""

from math import *
try:
    from numpy import *
except:
    try:
        from Numeric import *
    except:
        print "Could not import either Numeric nor numpy."
from copy import copy

#-----------------------------------------------------------------------------
def sinh(z):
    """
    Returns the hyperbolic sine of z
    """
    return (exp(z) - exp(-z)) / 2.0

def cosh(z):
    """
    Returns the hyperbolic cosine of z
    """
    return (exp(z) + exp(-z)) / 2.0

def sech(z):
    """
    Returns the hyperbolic secant or inverse cosh
    """
    return 1.0 / cosh(z)

def tanh(z):
    """
    Returns the hyperbolic tangent of z
    """
    return sinh(z) / cosh(z)

#-----------------------------------------------------------------------------

def find_beta(s, n):
    """
    Returns beta value from tanh(beta * x) which generates a suitable
    hyperbolic tangent curve for clustering nodes. 
    
    @param s: The user specified spacing for one end of the line
    @type s: float
    @param end: The end of the line for which the spacing is desired
    @type end: Integer (0=start of line, 1=end of line)
    @param n: number of nodes along interval
    @type n: Integer
    @return: This function returns the scaling factor beta which is applied
        to the function tanh( beta * x )
    """
    # We want n values in the (mathematical) range [0.0, 1.0]
    dx = 1.0/(n-1)
    x = array(range(0,n), Float) * dx
#     Commented 09-Dec-2005 - fixing line flip
    if 1:
        x1 = x[-2]; x2 = x[-1]
    else:
        x1 = x[0]; x2 = x[1]
    # evaluate s(beta) = tanh(beta*x2) - tanh(beta*x1) at both ends
    beta = arange(0, 10, 0.0001)
    sbeta_s0 = tanh(beta * x[1]) - tanh(beta * x[0])
    sbeta_s1 = tanh(beta * x[-1]) - tanh(beta * x[-2])
    # for the range of beta values, determine what the max possible spacing is
    sbeta_s0_max = sbeta_s0[argmax(sbeta_s0)]
    sbeta_s1_max = sbeta_s1[argmax(sbeta_s1)]
    # if max spacing is less than the user wanted, set as spacing to use
    if sbeta_s0_max < s and sbeta_s1_max < s:
        print "\nDesired spacing greater than both end max values.\n"
        if sbeta_s0_max > sbeta_s1_max:
            print "End 0 spacing largest. High resolution end=1.\n"
            s = sbeta_s0_max
            nominal_beta = beta[argmax(sbeta_s0)]
            flip = 1
        else:
            print "End 1 spacing largest. High resolution end=0.\n"
            s = sbeta_s1_max
            nominal_beta = beta[argmax(sbeta_s1)]
            flip = 0
    else:
        print "\nCan match desired spacing."
#--------------------------------------
#         Commented 09-Dec-2005 - see VERSION HISTORY
#         nominal_beta = beta[argsort(sbeta)[searchsorted(sort(sbeta), s)]]
#         s = sort(sbeta)[searchsorted(sort(sbeta), s)]
#         if abs(delta_s0 - s) > abs(delta_s1 - s):
#--------------------------------------
        # Checks if either max possible spacing is less than desired value
        # and if so, sets extremly large space for elimination.
        # Otherwise, sets nominal beta and space values for each end.
        if sbeta_s0_max < s:
            s0 = 1e10
            nominal_beta_s0 = None
        else:
            nominal_beta_s0 = beta[argsort(sbeta_s0)[searchsorted(sort(sbeta_s0), s)]]
            s0 = sort(sbeta_s0)[searchsorted(sort(sbeta_s0), s)]
        if sbeta_s1_max < s:
            s1 = 1e10
            nominal_beta_s1 = None
        else:
            nominal_beta_s1 = beta[argsort(sbeta_s1)[searchsorted(sort(sbeta_s1), s)]]
            s1 = sort(sbeta_s1)[searchsorted(sort(sbeta_s1), s)]

#         print "spacings ", s0, s1
#         print "betas ", nominal_beta_s0, nominal_beta_s1
        # Compare calculated spacing at each end and select closest to user
        # specified. In most cases, it will be spot on (within tolerance)
#         print "s = ", s
#         print "s0 = ", s0
#         print "s1 = ", s1
        if abs(s0 - s) < abs(s1 - s):
            print "End 0 spacing closest to desired. High resolution end=1."
            nominal_beta = nominal_beta_s0
            flip = 1
        else:
            print "End 1 spacing closest to desired. High resolution end=0."
            nominal_beta = nominal_beta_s1
            flip = 0
    # print "optimal beta=", nominal_beta
    return (nominal_beta, flip)

#-----------------------------------------------------------------------------

def smooth_curve(xp, yp, iters):
    """
    This function provides a means of smoothing a curve by assuming a linear
    gradient through node from its two neighbours.

    @param xp: Array of x values
    @type xp: Array
    @param yp: Array of y values
    @type yp: Array
    @param iters: Number of iterations over which to smooth the line
    @type iters: Integer
    """

    ystar = yp
    # iterate through max number of iterations and perform smoothing 
    # operations
    for i in range(0, iters):
        for j in range(1, len(yp)-1):
            ystar[j] = (yp[j+1] + yp[j-1]) / 2.0
        yp = ystar
    return ystar

#-----------------------------------------------------------------------------
def make_cluster(length, n, s0, s1, end0=None, end1=None, iters=10, midFrac=0.2):
    """
    This function provides the method for making the cluster array.

    @param length: The length of the block edge to be manipulated.
    @type length: Float
    @param n: Number of nodes to be distributed along the interval
    @type n: Integer
    @param s0: Spacing (dimensioned) desired at end 0 of the interval
    @type s0: Float
    @param s1: Spacing (dimensioned) desired at end 1 of the interval
    @type s1: Float
    @param iters: Number of iterations used for smoothing function 
    (Default: 5)
    @type iters: Integer
    @param midFrac: Fraction of nodes in middle of interval to be used for
    smoothing (Default: 0.2)
    @type midFrac: Float
    """
    s0 = s0 / length * 2.0
    s1 = s1 / length * 2.0
    # Split the collection into left and right pieces
    # keeping the same number of points overall.
    n0 = int(0.5 * n)
    n1 = n - n0

    # determine cluster parameters for both left and right halves of interval

    zeta0 = array(range(0,n0), Float)
    zeta0 = zeta0 / zeta0[-1]
    beta0, flip0 = find_beta(s0, n0)

    zeta1 = array(range(0,n1), Float)
    zeta1 = zeta1 / zeta1[-1]
    beta1, flip1 = find_beta(s1, n1)

    if abs(beta0-beta1) > 0.1:
        if beta0 > beta1:
            beta1 = beta0
        else:
            beta0 = beta1

    # perform clustering operations on left half of the interval
#     zeta0 = array(range(0,n0), Float)
#     zeta0 = zeta0 / zeta0[-1]
#     beta, flip = find_beta(s0, n0)
#     print "Beta = ", beta0
#     print "Length = ", length
#     print "n0 = ", n0
    etabar0 = tanh(zeta0 * beta0)
    etabar0 = etabar0 / etabar0[-1]
    # Perform flip operation if required 
    if 1: #flip == 0 or end0 == 1:
        # flip first half of interval and shift it so that it starts at zero
        left_piece = (etabar0 * -1.0) + 1.0
        left_piece = left_piece[::-1]
    else:
        left_piece = etabar0
    print 's0 achieved = ', 0.5 * (left_piece[1] - left_piece[0]) * length

    # Right-hand piece of original interval.
#     zeta1 = array(range(0,n1), Float)
#     zeta1 = zeta1 / zeta1[-1]
#     beta, flip = find_beta(s1, n1)
#     print "Beta = ", beta1
#     print "Length = ", length
#     print "n0 = ", n1
    etabar1 = tanh(zeta1 * beta1)
    etabar1 = etabar1 / etabar1[-1]
    # Perform flip operation if required 
    if 0:#flip == 1 or end1 == 0:
        # flip second half of interval
        right_piece = (etabar1 * -1.0) + 1.0
        right_piece = right_piece[::-1] + 1.0
    else:
        # offset the right_piece so that it starts where the left-piece finished.
        right_piece = etabar1 + 1.0
    print 's1 achieved = ', 0.5 * (right_piece[-1] - right_piece[-2]) * length
    # print "0. len(L)=", len(left_piece), " len(R)=", len(right_piece)
    #
    # combine two halves of the interval into one array
    etabar_final = concatenate( (left_piece, right_piece) )
    etabar_final /= etabar_final[-1]
    # print "1. len(etabar_final)=", len(etabar_final)

    # perform smoothing operation on middle fraction of points
    zeta_final = concatenate( (zeta0, zeta1+zeta0[-1]+1.0/n) ) / 2.0
    # print "1. len(zeta_final)=", len(zeta_final)
    midIdxs = \
        range(int(round(len(etabar_final)/2-midFrac/2*len(etabar_final))), \
        int(round(len(etabar_final)/2+midFrac/2*len(etabar_final))))
    etamid = smooth_curve(take(zeta_final, midIdxs),
                          take(etabar_final,midIdxs), iters)
    # replace middle fraction of smoothed points
    put(etabar_final, midIdxs, etamid)
    #
    # return an array of clustered points from 0.0 to 1.0
    # Fix rounding error.
    if etabar_final[0] < 0.0: etabar_final[0] = 0.0
    if etabar_final[-1] > 1.0: etabar_final[-1] = 1.0
    # print "2. len(etabar_final)=", len(etabar_final)
    return etabar_final


#-----------------------------------------------------------------------------
if __name__ == '__main__':
    from pylab import *
    length = 50.0e-3
    n = 100
    s0 = 1.0e-6
    s1 = 1.0e-3
    iters = 10
    midFrac = 0.6

    cluster = make_cluster(length, n, s0, s1, iters, midFrac)
    print "n=", n, "len(cluster)=", len(cluster)
    print "first=", cluster[0], "second=", cluster[1], \
          "diff=", cluster[1] - cluster[0]
    print "next-to-last=", cluster[-2], "last=", cluster[-1], \
          "diff=", cluster[-1] - cluster[-2]
    plot(cluster,'.')
    savefig('hyptan',dpi=150)
    show()
