# analytic_solution.py
# Python version of the analytic solution described in Appendix A of
# C.J. Roy, C.C. Nelson, T.M. Smith and C.C. Ober
# Verification of Euler/Navier-Stokes codes using the method
# of manufactured solutions.
# Int J for Numerical Methods in Fluids 2004; 44:599-620
#
# PJ, 28-May-2011
# It essentially Rowan's code with more and renamed variables
# to bring it closer to the original paper.
# PJ, 30-June-2012
# Scale the disturbance to reduce its magnitude away from the centre.
# RJG, 06-June-2014
# Re-worked completely to use sympy
#
# Jianyong Wang, 13-July-2016
# Add the analytic solutions for k-w turbulence model verification.
# Reference: C.J. Roy, E. Tendean, S.P. Veluri, and R. Rifki
#            E.A. Luke and S. Hebert 
# Verification of RANS turbulence models in Loci-CHEM using the method of manufactured solutions.
# 18th AIAA Computational Fluid Dynamics Conference, AIAA

from sympy import *
R_air = 287.0
# Read case no.
fp = open('case.txt', 'r');
case_str = fp.readline().strip()
case = int(case_str)
fp.close()
# constants
L = 1.0
if case == 1 or case == 3:
    # Supersonic flow
    rho0=1.0; rhox=0.15; rhoy=-0.1; rhoxy=0.0; arhox=1.0; arhoy=0.5; arhoxy=0.0;
    u0=800.0; ux=50.0; uy=-30.0; uxy=0.0; aux=1.5; auy=0.6; auxy=0.0;
    v0=800.0; vx=-75.0; vy=40.0; vxy=0.0; avx=0.5; avy=2.0/3; avxy=0.0;
    p0=1.0e5; px=0.2e5; py=0.5e5; pxy=0.0; apx=2.0; apy=1.0; apxy=0.0

if case == 2 or case == 4:
    # Subsonic(turbulent) flow
    rho0=1.0; rhox=0.15; rhoy=-0.1; rhoxy=0.08; arhox=0.75; arhoy=1.0; arhoxy=1.25;
    u0=70.0; ux=7.0; uy=-8.0; uxy=5.5; aux=1.5; auy=1.5; auxy=0.6;
    v0=90.0; vx=-5.0; vy=10.0; vxy=-11.0; avx=1.5; avy=1.0; avxy=0.9;
    p0=1.0e5; px=0.2e5; py=0.175e5; pxy=-0.25e5; apx=1.0; apy=1.25; apxy=0.75;
    tke0=780.0; tkex=160.0; tkey=-120.0; tkexy=80.0; atkex=0.65; atkey=0.7; atkexy=0.8;
    omega0=150.0; omegax=-30.0; omegay=22.5; omegaxy=40.0; a_omegax=0.75; a_omegay=0.875; a_omegaxy=0.6

x, y, rho, u, v, p, tke, omega = symbols('x y rho u v p tke omega')

rho = rho0 + rhox*cos(arhox*pi*x/L) + rhoy*sin(arhoy*pi*y/L) + \
   rhoxy*cos(arhoxy*pi*x*y/(L*L));
u =  u0 + ux*sin(aux*pi*x/L) + uy*cos(auy*pi*y/L) + uxy*cos(auxy*pi*x*y/(L*L));
v =  v0 + vx*sin(avx*pi*x/L) + vy*cos(avy*pi*y/L) + vxy*cos(avxy*pi*x*y/(L*L));
p =  p0 + px*cos(apx*pi*x/L) + py*sin(apy*pi*y/L) + pxy*sin(apxy*pi*x*y/(L*L));
tke =  tke0 + tkex*cos(atkex*pi*x/L) + tkey*sin(atkey*pi*y/L) + tkexy*cos(atkexy*pi*x*y/(L*L));
omega = omega0 + omegax*cos(a_omegax*pi*x/L) + omegay*sin(a_omegay*pi*y/L) + omegaxy*cos(a_omegaxy*pi*x*y/(L*L));

def ref_function(x1, y1, z1, t):
    inp = {x:x1, y:y1}
    rho1 = rho.subs(inp).evalf()
    p1 = p.subs(inp).evalf()
    T1 = p1/(rho1*R_air)
    u1 = u.subs(inp).evalf()
    v1 = v.subs(inp).evalf()
    tke1 = tke.subs(inp).evalf()
    omega1 = omega.subs(inp).evalf()
    return {"rho":rho1, "p":p1, "T":T1, "vel.x":u1, "vel.y":v1, "tke":tke1, "omega":omega1}

if __name__ == "__main__":
    pt = {x:0.5, y:0.5}
    print 'rho=', rho.subs(pt).evalf(), \
        'u=', u.subs(pt).evalf(), \
        'v=', v.subs(pt).evalf(), \
        'p=', p.subs(pt).evalf(), \
        'tke=', tke.subs(pt).evalf(), \
        'omega=', omega.subs(pt).evalf()
        
