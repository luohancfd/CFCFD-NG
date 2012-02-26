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
   * Oblique-shock relations.
"""

from math import *
from ..nm.secant_method import solve

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
    print "Done."
    return
