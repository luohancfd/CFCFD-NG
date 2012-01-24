#!/usr/bin/env python
#
# Author: Brendan O'Flaherty
# Date: 25 Aug 2009
#
# Writes an analytical solution to Sod's shock tube

from math import sqrt, pow
from numpy import sign, zeros
import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from cfpylib.nm.zero_finding import bisection_root

def p2p1_solver(p2p1, args):
    g1,p1,r1,g4,p4,r4 = args
    a1 = sqrt(g1*p1/r1)
    a4 = sqrt(g4*p4/r4)
    
    T2 = (g4 - 1.0)*(a1/a4)*(p2p1 - 1.0)/sqrt(2.0*g1*(2.0*g1 + (g1 + 1.0)*(p2p1 - 1.0)))
    return ((p4/p1) - p2p1*pow((1.0 - T2), -2.0*g4/(g4 - 1.0)))

def main():
    # Properties left of the diaphragm
    r4 = 1.0
    p4 = 1.0e5
    g4 = 1.4
    R4 = 288.2
    
    # Properties right of the diaphragm
    r1 = 1.0
    p1 = 1.0e4
    g1 = 1.4
    R1 = 288.2

    # Enter in at what time you want the profile.
    t = 1.0e-3

    # Enter in diaphragm location
    x0 = 0.5

    # How long do you want the shock tube to be?
    L = 1.0

    # How many points do you want displayed in the expansion wave?
    NX = 100
  
    # Calculations begin - no need to enter in values here
    T1 = p1/(R1*r1)
    T4 = p4/(R4*r4)

    u1 = 0
    u4 = 0

    a1 = sqrt(g1*p1/r1)
    a4 = sqrt(g4*p4/r4)
    
    args0=[g1,p1,r1,g4,p4,r4]

    if (x0 - a4*t < 0.0):
        print("Expansion wave reaches end of tube")
        print("Minimum diaphragm position must be: %e" % (a4*t))
        return 1
    if (L - (x0 + a1*t) < 0.0):
        print("Shock wave reaches end of tube")
        print("Minimum length must be: %e" % (x0 + a1*t))
        return 1
    if(sign(p2p1_solver(0, args0)) == sign(p2p1_solver(p4/p1, args0))):
        print("p2p1_solver(0) = %g, (p4/p1) = %g", 
              p2p1_solver(0, args0),p2p1_solver(p4/p1, args0))
        print("Please try again - signs are the same")
        return 1

    TOL=1e-9
    # These limits should always be sufficient to guarantee a root
    p2p1 = bisection_root(p2p1_solver, 0, p4/p1, TOL, args0)
    print("p2/p1 = %12.11e" % (p2p1))
    
    r2r1 = (1.0 + ((g1+1.0)/(g1-1.0))*(p2p1))/((g1+1.0)/(g1-1.0) + p2p1)
    W = a1*sqrt(((g1+1.0)/(2*g1))*(p2p1-1.0) + 1.0)
    up = W*(1 - (1/r2r1))

    print("r2/r1 = %12.11e, W = %12.11e, up = %12.11e" % (r2r1, W, up))

    r2 = r1*r2r1
    p2 = p1*p2p1
    T2 = p2/(R1*r2)
    u2 = up
    xs = x0+W*t # Location of shock
    xc = x0+up*t # Location of contact discontinuity
  
    p3p4 = p2p1/(p4/p1)
    r3r4 = pow(p3p4, 1.0/g4)
  
    r3 = r3r4*r4
    p3 = p3p4*p4
    a3 = sqrt(g4*p3/r3)
    u3 = up
    T3 = p3/(r3*R4)
  
    dx = ((t*(up-a3)) - (-a4*t))/NX
    rx = zeros([NX,2], float)
    px = zeros([NX,2], float)
    Tx = zeros([NX,2], float)
    ux = zeros([NX,2], float)

    rx[0][0] = x0 - a4*t
    px[0][0] = rx[0][0]; Tx[0][0] = rx[0][0]; ux[0][0] = rx[0][0];
    for i in range(1,NX):
        rx[i][0] = rx[i-1][0] + dx
        px[i][0] = rx[i][0]
        Tx[i][0] = rx[i][0]
        ux[i][0] = rx[i][0]
  
    for i in range(NX):
        u = (2.0/(g4+1.0))*(a4 + (rx[i][0]-x0)/t)
        rx[i][1] = r4*pow((1.0-((g4-1.0)/2.0)*(u/a4)), 2.0/(g4-1.0))
        ux[i][1] = u
        px[i][1] = p4*pow((1.0-((g4-1.0)/2.0)*(u/a4)), 2.0*g4/(g4-1.0))
        Tx[i][1] = px[i][1]/(R4*rx[i][1])
  
    fout = open("sod-analytical.dat","w")
    fout.write("# x, rho, p, T, u\n")
 
    # We go left to right
    # Region 4 to expansion fan
    fout.write("# region 4 start to expansion fan\n")
    fout.write("%e %e %e %e %e\n" % (0.0, r4, p4, T4, u4))
    fout.write("%e %e %e %e %e\n" % (x0-a4*t, r4, p4, T4, u4))
    
    # Expansion fan
    fout.write("# expansion fan\n")
    for i in range(NX):
        fout.write("%e %e %e %e %e\n" % (rx[i][0], rx[i][1], px[i][1], Tx[i][1], ux[i][1]))
    
    # Fan to contact surface
    fout.write("# fan to contact surface\n")
    fout.write("%e %e %e %e %e\n" % (t*(up-a3)+x0, r3, p3, T3, u3))
    fout.write("%e %e %e %e %e\n" % (xc, r3, p3, T3, u3))

    # Contact surface to shock
    fout.write("# contact surface to shock\n")
    fout.write("%e %e %e %e %e\n" % (xc, r2, p2, T2, u2))
    fout.write("%e %e %e %e %e\n" % (xs, r2, p2, T2, u2))

    # Shock to region 1's end
    fout.write("# shock to region 1 end\n")
    fout.write("%e %e %e %e %e\n" % (xs, r1, p1, T1, u1))
    fout.write("%e %e %e %e %e\n" % (L, r1, p1, T1 ,u1))

    fout.close() 

    return 0

if __name__=='__main__':
    main() 
