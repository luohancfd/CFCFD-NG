# Author: Rowan J. Gollan
# Place: The University of Queensland, Brisbane, Australia
# Date: 06-Jun-2014
#
# This script is used to generate the analytical source
# terms required to run the Method of Manufactured Solutions
# test case. The generated code is in Fortran95 and it can
# be converted to Lua with a separate script.
#
# 13-July-2016 (Jianyong Wang)
# The source term generations for k-w turbulent flow equations
# are implemented. Specific attentions have been paid in the
# max/min functions, and conditional statement.
# *** TURBULENCE EQUATIONS *****
# W.Y.K. Chan, P.A. Jacobs, J.P. Nap, D.J. Mee, R.M. Kirchhartz and S.J. Stennett
# The k-w turbulence model in Eilmer3: User guide and test cases.
# The University of Queensland, School of Mechanical & Mining Engineering, 
# Research Report Number 2010/01, pg 04-05.
# 
# This is an exercise in using sympy to generate the source
# terms. It is a transliteration of PJ's original work
# done using Maxima.

from sympy import *
from analytic_solution import *
import math

Rgas, g, Prandtl, PrT, Cv, Cp = symbols('Rgas g Prandtl PrT Cv Cp')
Rgas = 287.0
g = 1.4
Prandtl = 1.0
PrT = 0.89
Cv = Rgas/(g-1)
Cp = g*Cv

mu, k = symbols('mu k')
mu = 10.0
k = Cp*mu/Prandtl
if case == 1 or case == 3:
    mu = 0.0
    k = 0.0

# Thermodynamic behvaiour, equation of state and energy equation
e, T, et, ht = symbols('e T et ht')
e = p/rho/(g-1)
T = e/Cv
et = e + u*u/2 + v*v/2
ht = et + p/rho

# K_omega estimate of the turbulence viscosity
C_lim, beta_star, S_bar_squared, omega_t, mu_t = symbols('C_lim beta_star S_bar_squared omega_t mu_t')
C_lim = 0.875
beta_star = 0.09
S_bar_squared = diff(u, x)*diff(u, x) + diff(v, y)*diff(v, y) - 4.0/9.0*(diff(u, x) + diff(v, y))*(diff(u, x) + diff(v, y)) + \
                0.5*(diff(u, y) + diff(v, x))*(diff(u, y) + diff(v, x))
S_bar_squared = Max(0.0, S_bar_squared)
omega_t = Max(omega, C_lim*(2.0*S_bar_squared/beta_star)**0.5)
mu_t = rho*tke/omega_t

# Laminar and turbulent heat flux terms
qlx, qly, qtx, qty = symbols('qlx qly qtx qty')
qlx = -k*diff(T, x)
qly = -k*diff(T, y)
qtx = -Cp*mu_t*diff(T, x)/PrT
qty = -Cp*mu_t*diff(T, y)/PrT

# Laminar stress tensor
txx, tyy, txy = symbols('txx tyy txy')
txx = 2.0*mu*diff(u, x) - 2.0/3.0*mu*(diff(u, x) + diff(v, y))
tyy = 2.0*mu*diff(v, y) - 2.0/3.0*mu*(diff(u, x) + diff(v, y))
txy = mu*(diff(v, x) + diff(u, y))

# Turbulent stress tensor
tauxx, tauyy, tauxy = symbols('tauxx tauyy tauxy')
tauxx = 2.0*mu_t*diff(u, x) - 2.0/3.0*(mu_t*(diff(u, x) + diff(v, y)) + rho*tke)
tauyy = 2.0*mu_t*diff(v, y) - 2.0/3.0*(mu_t*(diff(u, x) + diff(v, y)) + rho*tke)
tauxy = mu_t*(diff(v, x) + diff(u, y))

# Closure coefficients for 2D cartesian flows
alpha, beta_0, beta, sigma, sigma_star, P_K, D_K, P_W, D_W, cross_diff, sigma_d, WWS, X_w, f_beta = symbols('alpha beta_0 beta sigma sigma_star P_K D_K P_W D_W cross_diff sigma_d WWS X_w f_beta')
alpha = 0.52
beta_0 = 0.0708
sigma = 0.5
sigma_star = 0.6
D_K = beta_star*rho*tke*omega
P_K = 1.3333*mu_t*(diff(u, x)*diff(u, x) - diff(u, x)*diff(v, y) + diff(v, y)*diff(v, y)) + \
   mu_t*(diff(u, y) + diff(v, x))*(diff(u, y) + diff(v, x)) - \
   0.66667*rho*tke*(diff(u, x) + diff(v, y))
P_K = Min(P_K, 25.0*D_K)
WWS = 0.0
cross_diff = diff(tke, x)*diff(omega, x) + diff(tke, y)*diff(omega, y)
sigma_d = Min(0.0, cross_diff)
P_W = alpha*omega/Max(0.1, tke)*P_K + sigma_d*rho/Max(1.0, omega)*cross_diff
X_w = abs(WWS/(beta_star*omega)**3.0)
f_beta = (1.0 + 85.0*X_w)/(1.0 + 100.0*X_w)
beta = beta_0*f_beta
D_W = beta*rho*omega*omega

# Wilcox k-w model equations for compressible flows in conservative form
t, fmass, fxmom, fymom, fe, ftke, fomega = symbols('t fmass fxmom fymom fe ftke fomega')
fmass = diff(rho, t) + diff(rho*u, x) + diff(rho*v, y)
fxmom = diff(rho*u, t) + diff(rho*u*u, x) + diff(rho*u*v, y) + diff(p, x) - diff(txx+tauxx, x) - diff(txy+tauxy, y)
fymom = diff(rho*v, t) + diff(rho*v*u, x) + diff(rho*v*v, y) + diff(p, y) - diff(txy+tauxy, x) - diff(tyy+tauyy, y)
fe = diff(rho*(et+tke), t) + diff(rho*u*(ht+tke), x) + diff(rho*v*(ht+tke), y) + diff(qlx+qtx, x) + diff(qly+qty, y) - \
     diff(u*(txx+tauxx)+v*(txy+tauxy), x) - diff(u*(txy+tauxy)+v*(tyy+tauyy), y) - \
     diff((mu+sigma_star*rho*tke/omega)*diff(tke, x), x) - diff((mu+sigma_star*rho*tke/omega)*diff(tke, y), y)
ftke = diff(rho*tke, t) + diff(rho*u*tke, x) + diff(rho*v*tke, y) - P_K + D_K - diff((mu+sigma_star*rho*tke/omega)*diff(tke, x), x) - \
       diff((mu+sigma_star*rho*tke/omega)*diff(tke, y), y)
fomega = diff(rho*omega, t) + diff(rho*u*omega, x) + diff(rho*v*omega, y) - P_W + D_W - diff((mu+sigma*rho*tke/omega)*diff(omega, x), x) - \
         diff((mu+sigma*rho*tke/omega)*diff(omega, y), y)

if __name__ == '__main__':
    import re
    from sympy.utilities.codegen import codegen
    print 'Generating manufactured source terms.'
    [(f_name, f_code), (h_name, f_header)] = codegen(
        [("fmass", fmass), ("fxmom", fxmom), ("fymom", fymom), ("fe", fe), ("ftke", ftke), ("fomega", fomega)],
        "F95", "test", header=False)
    # Convert F95 to Lua code
    # This is heavily borrowed PJ's script: f90_to_lua.py
    # First we'll do some replacements
    f_code = f_code.replace('**', '^')
    f_code = f_code.replace('d0', '')
    f_code = f_code.replace('Max', 'max')
    f_code = f_code.replace('Min(0.0', 'Min1(0.0')
    f_code = f_code.replace('Min(25.0', 'Min2(25.0')
    # Now we'll break into lines so that we can completely remove
    # some lines and tidy others
    lines = f_code.split('\n')
    lines[:] = [l.lstrip() for l in lines if ( not l.startswith('REAL*8') and
                                               not l.startswith('implicit') and
                                               not l.startswith('end')) ]
    # Now reassemble but collect the split lines into a large line
    buf = ""
    f_code = ""
    for i,l in enumerate(lines):
        if l.endswith('&'):
            buf = buf + l[:-2]
        else:
            if buf == "":
                f_code = f_code + l + '\n'
            else:
                f_code = f_code + buf + l + '\n'
                buf = ""

    fin = open('udf-source-template.lua', 'r')
    template_text = fin.read()
    fin.close()
    lua_text = template_text.replace('<insert-source-terms-here>',
                                     f_code)

    fout = open('udf-source.lua', 'w')
    fout.write(lua_text)
    fout.close()
    print 'Done converting to Lua.'

    



