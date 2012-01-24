#!/usr/bin/env python
## \file valliammai.py
## \brief Python script demonstrating grid clustering developed by Valliammai
## \author Adriaan Window
## \version January 2005
##
"""
This script generates a clustered node distribution as prescribed by the 
following paper:

Valliammai, V., Gogoi, A., Grid Quality Improvement using Multi-Block 
Smoothing, Computational Fluid Dynamics Journal, v. 10 no. 2, pp. 169-174, 
July 2001
"""

try:
    from numpy import *
except:
    try:
        from Numeric import *
    except:
        print "Could import neither numpy nor Numeric."



def smooth_curve(yp, iters):
    """
    This function provides a means of smoothing a curve by assuming a linear
    gradient through node from its two neighbours.

    @param yp: Array of y values
    @type yp: Array
    @param iters: Number of iterations over which to smooth the line
    @type iters: Integer
    """
    ystar = yp
    for i in range(0, iters):
        for j in range(1, len(yp)-1):
            ystar[j] = (yp[j+1] + yp[j-1]) / 2.0
        yp = ystar
    return ystar

def make_cluster(length, n, s0, s1, iters=10, midFrac=0.2):
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
    n0 = int(0.5 * n)
    n1 = n - n0

    zeta0 = array(range(0,n0), Float)
    zeta0 = zeta0 /zeta0[-1]
    alpha = find_alpha(s0, n0)
    etabar0 = valliammai(alpha, zeta0)
    etabar0 = etabar0 / etabar0[-1]
    left_piece = etabar0

    zeta1 = array(range(0,n1), Float)
    zeta1 = zeta1 /zeta1[-1]
    alpha = find_alpha(s1, n1)
    etabar1 = valliammai(alpha, zeta1)
    etabar1 = etabar1 / etabar1[-1]
    right_piece = etabar1
    right_piece = (etabar1 * -1.0) + 1.0
    right_piece = right_piece[::-1] + 1.0

    # join left and right halves
    etafinal = concatenate( (left_piece, right_piece) )
    etafinal /= etafinal[-1]
    if etafinal[0] < 0.0: etafinal[0] = 0.0
    if etafinal[-1] > 1.0: etafinal[-1] = 1.0

    # smooth middle of curve
    midIdxs = \
        range(int(round(len(etafinal)/2-midFrac/2*len(etafinal))), \
        int(round(len(etafinal)/2+midFrac/2*len(etafinal))))
    etamid = smooth_curve(take(etafinal, midIdxs), iters)
    put(etafinal, midIdxs, etamid)
    s0 = 0.5 * (left_piece[1] - left_piece[0]) * length
    s1 = 0.5 * (right_piece[-1] - right_piece[-2]) * length
    return etafinal, s0, s1


def find_alpha(s, n):
    """
    Returns alpha value from Valliammai function which generates a suitable
    curve for clustering nodes. 
    
    @param s: The user specified spacing for one end of the line
    @type s: float
    @param end: The end of the line for which the spacing is desired
    @type end: Integer (0=start of line, 1=end of line)
    @param n: number of nodes along interval
    @type n: Integer
    @return: This function returns the scaling factor alpha which is applied
        to the function tanh( beta * x )
    """
    dx = 1.0/(n-1)
    x = array(range(0,n), Float) * dx

    alpha_ary = arange(1, 10, 0.0001)
    salpha = arange(1, 10, 0.0001)
    salpha = valliammai(alpha_ary, x[1]) - valliammai(alpha_ary, x[0])
    pyalpha = open('pyalpha.txt','w');
    for i in range(0,len(salpha)):
        pyalpha.write(str(alpha_ary[i])+" "+str(salpha[i])+"\n");
    pyalpha.close()
    salpha_max = salpha[argmax(salpha)]
    if salpha_max < s:
        s = salpha_max
        nominal_alpha = alpha_ary[argmax(salpha)]
    else:
        nominal_alpha = alpha_ary[argsort(salpha)[searchsorted(sort(salpha), s)]]
    print nominal_alpha
    return nominal_alpha

def valliammai(alpha, t):
    """
    Returns an evaluation of the valliammai function.
    """
    return (exp(alpha * t) - 1.0) / (exp(alpha) - 1.0)

if __name__ == "__main__":
    from string import *
    length = 0.01
    n = 100
    s0 = 1.0e-5
    s1 = 1.0e-5
    iters = 10
    midFrac = 0.2

    # simple cluster plot test
#    from pylab import *
    cluster  = make_cluster(length, n, s0, s1, iters=10, midFrac=0.2)
#    plot(cluster,'.')
#    show()

    # variance test
#     import Gnuplot
#     length = arange(0.01, 0.1, 0.01)
#     pts = range(50,150)
#     s0ary = zeros((len(length),len(pts)), Float)
#     s1ary = zeros((len(length),len(pts)), Float)
#     s0file = open('s0ary.txt','w')
#     s1file = open('s1ary.txt','w')
#     for l in range(0, len(length)):
#         print l
#         for n in range(0, len(pts)):
#             cluster, s0out, s1out = make_cluster(length[l], pts[n], s0, s1, iters, midFrac)
#             s0ary[l][n] = s0out
#             s1ary[l][n] = s1out
#             s0file.write(str(s0out)+' ')
#             s1file.write(str(s1out)+' ')
#         s0file.write('\n')
#         s1file.write('\n')
#     s1file.close()
#     s0file.close()
#     g=Gnuplot.Gnuplot()
#     g('set data style lines')
#     g('set xlabel "Interval length"')
#     g('set ylabel "Node count"')
#     g('set zlabel "Achieved spacing"')
#     g('set title "Spurious instances of cluster spacing"')
#     g.splot(Gnuplot.GridData(s0ary, length, pts, binary=0))
# 
#     raw_input('Press a key...\n')
#     g.hardcopy('hyptan_test.ps', enhanced=1, color=1)
