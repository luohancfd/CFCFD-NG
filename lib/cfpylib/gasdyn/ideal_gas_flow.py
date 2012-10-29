"""
ideal_gas_flow.py: One-dimensional steady flow of an ideal gas.

.. Author:
   PA Jacobs
   Centre for Hypersonics, School of Engineering
   The University of Queensland

.. Versions:
   1.1 30-Sep-94: Xplore version
   2.0 16-May-04: Python equivalent adapted from the Xplore version.
   27-Feb-2012: use relative import in cfpylib

Contents:

* One-dimensional flows:

   * Isentropic flow relations.
     State zero (0) refers to the stagnation condition.
     State star is the sonic (throat) condition.
   * 1D (Normal) Shock Relations
     State 1 is before the shock and state 2 after the shock.
     Velocities are in a shock-stationary frame.
   * 1-D flow with heat addition (Rayleigh-line)
     State star is the (hypothetical) sonic condition.

* Two-dimensional flows:

   * Prandtl-Meyer functions
   * Oblique-shock relations
   * Taylor-Maccoll conical flow
"""

from math import *
import numpy
from ..nm.secant_method import solve
from ..nm.zero_solvers import secant

# ---------------------------------------------------------------
# Isentropic flow

def A_Astar(M, g=1.4):
    "Returns area ratio A/Astar for a specified Mach number."
    t1 = (g + 1.0) / (g - 1.0)
    m2 = M**2
    t2 = 1.0 / m2 * (2.0 / (g + 1.0) * (1.0 + (g - 1.0) * 0.5 * m2))**t1
    t2 = sqrt(t2)
    return t2

def T0_T(M, g=1.4):
    "Returns temperature ratio T0/T for a specified Mach number."
    return 1.0 + (g - 1.0) * 0.5 * M**2

def p0_p(M, g=1.4):
    "Returns pressure ratio p0/p for a specified Mach number."
    return (T0_T(M, g))**( g / (g - 1.0) )

def r0_r(M, g=1.4):
    "Returns density ratio r0/r for a specified Mach number."
    return (T0_T(M, g))**(1.0 / (g - 1.0))

# -----------------------------------------------------------------
# 1-D normal shock relations.

def m2_shock(M1, g=1.4):
    "Returns the Mach number M2 after a normal shock."
    numer = 1.0 + (g - 1.0) * 0.5 * M1**2
    denom = g * M1**2 - (g - 1.0) * 0.5
    return sqrt(numer / denom)

def r2_r1(M1, g=1.4):
    "Returns the density ratio r2/r1 across a normal shock."
    numer = (g + 1.0) * M1**2
    denom = 2.0 + (g - 1.0) *M1**2
    return numer / denom

def u2_u1(M1, g=1.4):
    "Returns the velocity ratio u2/u1 across a normal shock."
    return 1 / r2_r1(M1, g)

def p2_p1(M1, g=1.4):
    "Returns the pressure ratio p2/p1 across a normal shock."
    return 1.0 + 2.0 * g / (g + 1.0) * (M1**2 - 1.0)

def T2_T1(M1, g=1.4):
    "Returns the temperature ratio T2/T1 across a normal shock."
    return  p2_p1(M1, g) / r2_r1(M1, g)

def p02_p01(M1, g=1.4):
    "Returns the ratio of stagnation pressures p02/p01 across the shock."
    t1 = (g + 1.0) / (2.0 * g * M1**2 - (g - 1.0)) 
    t2 = (g + 1.0) * M1**2 / (2.0 + (g - 1.0) * M1**2)
    return t1**(1.0/(g-1.0)) * t2**(g/(g-1.0))

def DS_Cv(M1, g=1.4):
    "Returns the entropy change ds across the shock normalized by Cv."
    t1 = p2_p1(M1, g)
    t2 = r2_r1(M1, g)
    return log(t1 * t2**g)

def pitot_p(p1, M1, g=1.4):
    "Returns the pitot pressure for a specified Mach number. Will shock the gas if required."
    if M1 > 1.0:
        p2 = p2_p1(M1,g)*p1
        M2 = m2_shock(M1, g)
        return p0_p(M2, g)*p2
    else:
        return p0_p(M1, g)*p1


# -----------------------------------------------------------------
# 1-D flow with heat addition (Rayleigh-line)

def T0_T0star(M, g=1.4):
    "Returns the total temperature ratio T0/T0star for a given Mach number."
    term1 = (g + 1.0) * M**2
    term2 = (1.0 + g * M**2)**2
    term3 = 2.0 + (g - 1.0) * M**2
    return term1 / term2 * term3

def M_Rayleigh(T0T0star, g=1.4):
    "Computes M from Total Temperature ratio for Rayleigh-line flow."
    # supersonic flow is assumed for the initial guesses
    def f_to_solve(m): return T0_T0star(m, g) - T0T0star
    return solve(f_to_solve, 2.5, 2.4)

def T_Tstar(M, g=1.4):
    "Returns static temperature ratio T/Tstar for Rayleigh-line flow. "
    return M**2 * ( (1.0 + g) / (1.0 + g * M**2) )**2

def p_pstar(M, g=1.4):
    " Returns static pressure ratio p/pstar for Rayleigh-line flow."
    return (1.0 + g) / (1.0 + g * M**2)

def r_rstar(M, g=1.4):
    "Returns density ratio r/rstar for Rayleigh-line flow."
    return 1.0 / M**2 / (1.0 + g) * (1.0 + g * M**2)

def p0_p0star(M, g=1.4):
    "Returns stagnation pressure ratio p0/p0star for Rayleigh-line flow."
    term1 = (2.0 + (g - 1.0) * M**2) / (g + 1.0)
    term2 = g / (g - 1.0)
    return (1.0 + g) / (1.0 + g * M**2) * term1**term2

# -----------------------------------------------------------------
# Prandtl-Meyer functions

def deg_to_rad(d): return d / 180.0 * pi
def rad_to_deg(r): return r * 180.0 / pi

def PM1(M, g=1.4):
    "Returns the Prandtl-Meyer value (in radians) for a given Mach number."
    if M > 1.0:
        t1 = M**2 - 1.0
        t2 = sqrt((g - 1.0) / (g + 1.0) * t1)
        t3 = sqrt(t1)
        t4 = sqrt((g + 1.0) / (g - 1.0))
        nu = t4 * atan(t2) - atan(t3)
    else:
        nu = 0.0
    return nu

def PM2(nu, g=1.4):
    "Returns the Mach number from the Prandtl-Meyer value (in radians)."
    # Solve the equation PM1(m, g) - nu = 0, assuming supersonic flow.
    def f_to_solve(m): return PM1(m, g) - nu
    return solve(f_to_solve, 2.0, 2.1)

# -----------------------------------------------------------------
# Oblique shock relations
# beta is shock angle wrt on-coming stream direction (in radians)
# theta is flow deflection wrt on-coming stream (in radians)

def beta_obl(M1, theta, g=1.4):
    """Returns the oblique shock wave angle
    given the upstream Mach number and the deflection angle."""
    b1 = asin(1.0/M1); b2 = b1 * 1.05
    def f_to_solve(beta): return theta_obl(M1, beta, g) - theta
    return solve(f_to_solve, b1, b2)

def theta_obl(M1, beta, g=1.4):
    "Compute the deflection angle given the shock wave angle."
    m1sb = M1 * sin(beta)
    t1 = 2.0 / tan(beta) * (m1sb**2 - 1.0) 
    t2 = M1**2 * (g + cos(2.0 * beta)) + 2.0
    theta = atan(t1/t2)
    return theta

def M2_obl(M1, beta, theta, g=1.4):
    "Returns the Mach number after an oblique shock."
    m1sb = M1 * sin(beta)
    numer = 1.0 + (g - 1.0) * 0.5 * m1sb**2
    denom = g * m1sb**2 - (g - 1.0) * 0.5
    m2 = sqrt(numer / denom / (sin(beta - theta))**2 )
    return m2

def r2_r1_obl(M1, beta, g=1.4):
    "Compute the density ratio r2/r1 across an oblique shock"
    m1sb = M1 * sin(beta)
    numer = (g + 1.0) * m1sb**2
    denom = 2.0 + (g - 1.0) * m1sb**2
    return numer / denom

def u2_u1_obl(M1, beta, g=1.4):
    "Compute the flow-speed ratio u2/u1 across an oblique shock."
    return sqrt((sin(beta) / r2_r1_obl(M1, beta, g))**2 + (cos(beta))**2)

def p2_p1_obl(M1, beta, g=1.4):
    "Compute the pressure ratio p2/p1 across an oblique shock."
    m1sb = M1 * sin(beta)
    return 1.0 + 2.0 * g / (g + 1.0) * (m1sb**2 - 1.0)

def T2_T1_obl(M1, beta, g=1.4):
    "Compute the temperature ratio T2/T1 across a normal shock."
    return p2_p1_obl(M1, beta, g) / r2_r1_obl(M1, beta, g)

def p02_p01_obl(M1, beta, g=1.4):
    "Compute the ratio of stagnation pressures p02/p01 across an oblique shock."
    m1sb = M1 * sin(beta)
    t1 = (g + 1.0) / (2.0 * g * m1sb**2 - (g - 1.0)) 
    t2 = (g + 1.0) * m1sb**2 / (2.0 + (g - 1.0) * m1sb**2)
    return t1**(1.0/(g-1.0)) * t2**(g/(g-1.0))

#------------------------------------------------------------------------
# Taylor-Maccoll cone flow.

def taylor_maccoll_odes(z, theta, g=1.4):
    """
    The ODEs from the Taylor-Maccoll formulation.

    See PJ's workbook for Feb 2012 for details.
    We've packaged them formally so that we might one day use
    a more sophisticated ODE integrator requiring fewer steps.
    """
    rho, V_r, V_theta, h, p = z
    # Assemble linear system for determining the derivatives wrt theta.
    A = numpy.zeros((5,5), float)
    b = numpy.zeros((5,), float)
    A[0,0] = V_theta; A[0,2] = rho; b[0] = -2.0*rho*V_r - rho*V_theta/tan(theta)
    A[1,1] = 1.0; b[1] = V_theta
    A[2,1] = rho*V_r; A[2,2] = rho*V_theta; A[2,4] = 1.0
    A[3,1] = V_r; A[3,2] = V_theta; A[3,3] = 1.0
    A[4,0] = h*(g-1)/g; A[4,3] = rho*(g-1)/g; A[4,4] = -1.0
    dzdtheta = numpy.linalg.solve(A,b)
    return dzdtheta

def theta_cone(V1, p1, T1, beta, R=287.1, g=1.4):
    """
    Compute the cone-surface angle and conditions given the shock wave angle.

    :param V1: speed of gas into shock
    :param p1: free-stream pressure
    :param T1: free-stream static temperature
    :param beta: shock wave angle wrt stream direction (in radians)
    :param R: gas constant
    :param g: ratio of specific heats
    :returns: tuple of theta_c, V_c, p_c, T_c:
        theta_c is stream deflection angle in radians
        V_c is the cone-surface speed of gas in m/s
        p_c is the cone-surface pressure
        T_c is the cone-surface static temperature

    The computation starts with the oblique-shock jump and then integrates
    across theta until V_theta goes through zero.
    The cone surface corresponds to V_theta == 0.

    This ideal-gas version adapted from the cea2_gas_flow version, 08-Mar-2012.

    24-Jun-2012 : RJG added checks to catch the limiting case when beta < mu
                : and a linear interpolation when beta is only slightly larger
                : than mu (1% larger)
    """
    # When beta is only this fraction larger than mu,
    # we'll apply a linear interpolation
    LINEAR_INTERP_SWITCH = 1.01
    # Free-stream properties and gas model.
    a1 = sqrt(g*R*T1)
    M1 = V1 / a1
    C_p = R * g / (g-1)
    h1 = C_p * T1
    rho1 = p1 / (R * T1)
    # Test beta in relation to the Mach angle, mu
    mu = asin(1.0/M1)
    beta2 = LINEAR_INTERP_SWITCH*mu
    #print "beta= ", beta, "mu= ", mu, " beta2= ", beta2
    if beta <= mu:
        # An infinitely weak shock angle
        return 0.0, V1, p1, T1
    if beta < beta2:
        # It is difficult to integrate between the shock and cone body
        # when the shock angle is only slightly larger than the Mach
        # angle. In this instance, find the value at LINEAR_INTER_SWITCH*mu
        # and linearly interpolate to find the value at beta
        (theta2, V2, p2, T2) = theta_cone(V1, p1, T1, beta2, R, g)
        frac = (beta - mu)/(beta2 - mu)
        theta_c = frac*theta2
        V = (1.0 - frac)*V1 + frac*V2
        p = (1.0 - frac)*p1 + frac*p2
        T = (1.0 - frac)*T1 + frac*T2
        return theta_c, V, p, T
    #
    # Start at the point just downstream the oblique shock.
    theta_s = theta_obl(M1, beta, g)
    M2 = M2_obl(M1, beta, theta_s, g)
    assert M2 > 1.0
    rho2 = rho1 * r2_r1_obl(M1, beta, g)
    V2 = V1 * u2_u1_obl(M1, beta, g)
    p2 = p1 * p2_p1_obl(M1, beta, g)
    T2 = T1 * T2_T1_obl(M1, beta, g)
    h2 = T2 * C_p
    #
    # Initial conditions for Taylor-Maccoll integration.
    dtheta = -0.05 * pi / 180.0  # fraction-of-a-degree steps
    theta = beta
    V_r = V2 * cos(beta - theta_s)
    V_theta = -V2 * sin(beta - theta_s)
    # For integrating across the shock layer, the state vector is:
    z = numpy.array([rho2, V_r, V_theta, h2, p2])
    while V_theta < 0.0:
        # Keep a copy for linear interpolation at the end.
        z_old = z.copy(); theta_old = theta
        # Do the update using a low-order method (Euler) for the moment.
        dzdtheta = taylor_maccoll_odes(z, theta, g)
        z += dtheta * dzdtheta; theta += dtheta
        rho, V_r, V_theta, h, p = z
        if False: print "DEBUG theta=", theta, "V_r=", V_r, "V_theta=", V_theta
    # At this point, V_theta should have crossed zero so
    # we can linearly-interpolate the cone-surface conditions.
    V_theta_old = z_old[2]
    frac = (0.0 - V_theta_old)/(V_theta - V_theta_old)
    z_c = z_old*(1.0-frac) + z*frac
    theta_c = theta_old*(1.0-frac) + theta*frac
    # At the cone surface...
    rho, V_r, V_theta, h, p = z_c
    T = h / C_p
    assert abs(V_theta) < 1.0e-6
    #
    return theta_c, V_r, p, T

def beta_cone(V1, p1, T1, theta, R=287.1, g=1.4):
    """
    Compute the conical shock wave angle given the cone-surface deflection angle.

    :param V1: speed of gas into shock
    :param p1: free-stream pressure
    :param T1: free-stream static temperature
    :param theta: stream deflection angle (in radians)
    :param R: gas constant
    :param g: ratio of specific heats
    :returns: shock wave angle wrt incoming stream direction (in radians)

    This ideal-gas version adapted from the cea2_gas_flow version, 08-Mar-2012.
    """
    # Free-stream properties and gas model.
    a1 = sqrt(g*R*T1)
    M1 = V1 / a1
    C_p = R * g / (g-1)
    h1 = C_p * T1
    rho1 = p1 / (R * T1)
    # Initial guess
    M1 = V1 / a1
    b1 = asin(1.0 / M1) * 1.01 # to be stronger than a Mach wave
    b2 = b1 * 1.05
    def error_in_theta(beta_guess):
        theta_guess, V_c, p_c, T_c = theta_cone(V1, p1, T1, beta_guess, R, g)
        return theta_guess - theta
    return secant(error_in_theta, b1, b2, tol=1.0e-4, limits=[asin(1.0/M1), pi/2.0])

def beta_cone2(M1, theta, R=287.1, g=1.4):
    """
    Compute the conical shock wave angle given the cone-surface deflection angle and free stream Mach number.

    :param M1: free stream Mach number
    :param theta: stream deflection angle (in radians)
    :param R: gas constant
    :param g: ratio of specific heats
    :returns: shock wave angle wrt incoming stream direction (in radians)

    This version basically delegates work to beta_cone().
    """
    # Compute free stream velocity assuming unit value temperature
    T1 = 1.0
    a1 = sqrt(g*R*T1)
    V1 = M1*a1
    # Set free stream pressure to unit value
    p1 = 1.0
    # Now ready to call beta_cone()
    return beta_cone(V1, p1, T1, theta, R, g)

# -----------------------------------------------------------------

def demo():
    print "Begin test of isentropic flow ratios..."
    M = 2.0
    print "Computed: M=%g: A/Astar=%g, T0/T=%g, p0/p=%g, r0/r=%g" % \
          (M, A_Astar(M), T0_T(M), p0_p(M), r0_r(M))
    print "Expected: M=2, A/Astar=1.687, T0/T=1.80, p0/p=7.824, r0/r=4.347"
    print ""
    print "Normal shock jump..."
    print "Computed: M=%g: M2=%g, T2/T1=%g, p2/p1=%g, r2/r1=%g" % \
          (M, m2_shock(M), T2_T1(M), p2_p1(M), r2_r1(M))
    print "Expected: M1=2, M2=0.5774, T2/T1=1.687, p2/p1=4.50, r2/r1=2.667"
    print ""
    print "Rayleigh-line flow..."
    print "Computed: M=%g: T0/Tstar=%g, T/Tstar=%g, p/pstar=%g, r/rstar=%g" % \
          (M, T0_T0star(M), T_Tstar(M), p_pstar(M), r_rstar(M))
    print "Expected: M=2, T0/T0star=0.7934, T/Tstar=0.5289, p/pstar=0.3636, r/rstar=0.6875"
    print "Inverse calculation: T0/T0star=%g --> M=%g" % \
          (T0_T0star(M), M_Rayleigh(T0_T0star(M)))
    print ""
    print "Prandtl-Meyer function..."
    print "Computed: M=%g --> nu=%g; Inverse: M=%g <-- nu=%g" % \
          (M, PM1(M), PM2(1.1481), 1.1481)
    print "Expected: M=2 --> nu=0.4604; Inverse: M=4 <-- nu=1.1481"
    print ""
    print "Oblique shock relations may not quite match (data is from chart)..."
    beta = deg_to_rad(44.0); theta = deg_to_rad(14.0); # from chart, M=2
    print "Computed: M1=%g, theta(beta=%g)=%g, beta(theta=%g)=%g" % \
          (M, beta, theta_obl(M, beta), theta, beta_obl(M, theta))
    print "Conditions behind shock:"
    print "M2=%g, expected 1.482 (from chart, 14 degree deflection)" % \
          M2_obl(M, beta, theta)
    print "Computed: T2/T1=%g, p2/p1=%g, r2/r1=%g" % \
          (T2_T1_obl(M, beta), p2_p1_obl(M, beta), r2_r1_obl(M, beta))
    print "Expected: T2/T1=1.249, p2/p1=2.088, r2/r1=1.673 (arox. table M=1.390)"
    print "u2/u1=%g, p02/p01=%g" % \
          (u2_u1_obl(M, beta), p02_p01_obl(M, beta))
    print "Expected: u2/u1=0.8304=sin(B)/sin(B-d)*r1/r2"
    print ""
    M1 = 1.5; p1 = 100.0e3; T1 = 300.0; R = 287.1; g = 1.4; rho1 = p1/(R*T1)
    print "Taylor-Maccoll cone flow demo with M1=%g" % M1
    print "for M1=1.5, beta=49deg, expect theta=20deg from NACA1135."
    a1 = sqrt(1.4*287*T1)
    V1 = M1 * a1
    beta = 49.0 * pi/180
    theta_c, V_c, p_c, T_c = theta_cone(V1, p1, T1, beta)
    print "theta_c(deg)=", theta_c*180.0/pi, "expected 20deg, surface speed V_c=", V_c
    print "surface pressure coefficient=", (p_c - p1)/(0.5*rho1*V1*V1), "expected 0.385"
    print "p_c: %g, T_c: %g" % (p_c, T_c)
    print ""
    print "Conical shock from cone with half-angle 20deg in M1=", M1
    beta = beta_cone(V1, p1, T1, 20.0*pi/180)
    print "sigma(deg)=", beta*180/pi, "expected 49deg"
    print "Repeat above test, but call beta_cone2()"
    beta = beta_cone2(M1, 20.0*pi/180)
    print "sigma(deg)=", beta*180/pi, "expected 49deg"
    #
    print "Done."
    return
