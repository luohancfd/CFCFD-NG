#! /usr/bin/env python
# self_similar_bl.py
"""
Set up a boundary-layer flow profile.

The calculation is based on the self-similar solution for a flat plate,
as described in section 6.5 of J.D. Anderson's text
Hypersonic and High Temperature Gas Dynamics.

\author PJ
\version 06-Mar-2011
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
import math, numpy
import cfpylib.nm.ode as ode
import cfpylib.gasdyn.sutherland as sutherland
import cfpylib.nm.nelmin as nelmin

print "self_similar_bl: begin..."

print "Gas constants for ideal air"
R = 287.1  # J/kg.K
gamma = 1.4
C_p = gamma/(gamma-1.0) * R
Pr = 0.71  # to match Schetz' boundary-layer code
print "R=", R, "gamma=", gamma, "C_p=", C_p, "Pr=", Pr

print "Conditions just outside of boundary layer to match Schetz's worked example 5-1"
p_e = 1.013e3    # Pa
u_e = 1390.0   # m/s
T_e = 300.0    # degrees K
h_e = C_p * T_e
T_wall = 300.0
h_wall = C_p * T_wall
print "p_e=", p_e, "u_e=", u_e, "T_e=", T_e, "h_e=", h_e
print "Condition at wall"
print "T_wall=", T_wall, "h_wall=", h_wall

rho_e = p_e / (R * T_e) # ideal equation of state
mu_e = sutherland.mu(T_e,'Air')
k_e = mu_e * C_p / Pr
print "rho_e=", rho_e, "mu_e=", mu_e, "k_e=", k_e

# Choose a position along the plate
x = 0.1  # metres
xi = rho_e * u_e * mu_e * x
print "x=", x, "xi=", xi

def C(g):
    """
    Ratio of density.viscosity product at points in the boundary layer.
    """
    T = g * h_e / C_p
    T = max(T, 100.0) # to avoid difficult values
    rho = p_e / (R * T)
    mu = sutherland.mu(T,'Air')
    return rho*mu/(rho_e*mu_e)

def Cd(g, gd):
    """
    Finite-difference estimate of dC/deta.
    """
    deltag = 0.01 * g  # something not too big, not too small
    C0 = C(g)
    C1 = C(g+deltag)
    dCdg = (C1 - C0)/deltag
    return dCdg * gd

def odes(eta, z, n):
    """
    Functions defining the differential equations from Anderson's text.

    Elements of state vector z:
    f    stream function
    fd   normalized velocity df/deta = u/u_e
    fdd  d(fd)/deta
    g    normalized enthalpy h/h_e
    gd   dg/deta
    y    spatial coordinate through the boundary layer
    """
    assert len(z) == 6
    f, fd, fdd, g, gd, y = z
    fddd = 1.0/C(g) * (-f*fdd - Cd(g,gd)*fdd)
    gdd = Pr/C(g) * (-gd*(Cd(g,gd)/Pr+f) - C(g)*u_e*u_e/h_e*fdd*fdd)
    yd = math.sqrt(2*xi)/u_e * h_e/p_e * (gamma-1.0)/gamma * g
    dzdeta = numpy.array([fd, fdd, fddd, gd, gdd, yd])
    return dzdeta

def integrate_through_bl(fdd, gd):
    """
    Start at wall and integrate the differential equations through the BL.

    Input: fdd, gd: guesses for these elements
    Returns: eta and state vector values through the boundary layer (and beyond)
    """
    # Starting values at the wall.
    f = 0.0
    fd = 0.0
    g = h_wall/h_e
    y = 0.0
    z0 = numpy.array([f, fd, fdd, g, gd, y])
    # Integrate through the boundary layer to a large value of eta.
    eta = numpy.linspace(0.0, 5.0, 500)
    deta = eta[1] - eta[0]
    z_list = [z0]
    for i in range(1,len(eta)):
        eta1, z1, err = ode.rkf45_step(eta[i-1], deta, odes, 6, z_list[i-1])
        z_list.append(z1)
    return eta, z_list

def objective(params):
    """
    Evaluate the guessed values for fdd and gd by returning a measure
    of the error at the outer edge of the boundary layer.
    """
    assert len(params) == 2
    fdd, gd = params
    eta, z = integrate_through_bl(fdd, gd)
    f, fd, fdd, g, gd, y = z[-1]
    penalty_value = abs(fd - 1.0) + abs(g - 1.0)
    return penalty_value

if 0:
    # Determine the best initial values of fdd and gd.
    params, obj, flag, nfe, nrestart = nelmin.minimize(objective, [0.5, 1.0])
    print "params=", params, "obj=", obj, "nfe=", nfe, "nrestart=", nrestart
    # Integrate this particular boundary layer.
    fdd, gd = params
else:
    # I happened to prepare this solution a little earlier...
    # with obj= 3.22128965724e-08 nfe= 173 nrestart= 0
    fdd, gd = 0.44547389521432856, 1.0632176798422215
eta, z = integrate_through_bl(fdd, gd)

# Write the data for plotting (with gnuplot, maybe).
fp = open('profile.data', 'w')
fp.write("# eta f fd fdd g gd y p T rho u\n")
for i in range(len(eta)):
    f, fd, fdd, g, gd, y = z[i]
    h = g * h_e; T = h / C_p; rho = p_e / (R * T); u = fd * u_e
    fp.write("%f %f %f %f %f %f %f %f %f %f %f\n" %
             (eta[i], f, fd, fdd, g, gd, y, p_e, T, rho, u))
fp.close()
# A Lua table
fp = open('profile.lua', 'w')
fp.write("-- profile.lua: machine generated.\n")
fp.write("y={}; p={}; T={}; u={}\n")
for i in range(len(eta)):
    f, fd, fdd, g, gd, y = z[i]
    h = g * h_e; T = h / C_p; rho = p_e / (R * T); u = fd * u_e
    fp.write("y[%d]=%f; p[%d]=%f; T[%d]=%f; u[%d]=%f\n" %
             (i, 0.44-y, i, p_e, i, T, i, u))
fp.close()
# Python lists
fp = open('profile.py', 'w')
fp.write("# profile.py: machine generated.\n");
fp.write("y=[]; p=[]; T=[]; u=[]\n")
for i in range(len(eta)):
    f, fd, fdd, g, gd, y = z[i]
    h = g * h_e; T = h / C_p; rho = p_e / (R * T); u = fd * u_e
    fp.write("y.append(%f); p.append(%f); T.append(%f); u.append(%f)\n" % (0.44-y, p_e, T, u))
fp.close()

print "Done."
