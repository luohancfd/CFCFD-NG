#!/usr/bin/env python

from math import sqrt
from scipy.special import erfc
from libprep3 import *

#
# Gas setup
#

gmodel = create_gas_model("gas-model.lua")
Q = Gas_data(gmodel)

D = 9.88072e-06
p = 100.0e3
T = 273.2
Q.p = p
Q.T[0] = T
Q.massf[0] = 1.0
Q.massf[1] = 0.0
gmodel.eval_thermo_state_pT(Q)
rho_A = Q.rho

Q.massf[0] = 0.0
Q.massf[1] = 1.0
gmodel.eval_thermo_state_pT(Q)
rho_B = Q.rho

#
# geometry
#

xmin = -2.0e-5
xmax = 2.0e-5
dx = 1.0e-6

def calc_rho(r, t):
    arg = r / (2.0*sqrt(D*t))
    # erfc - error function complement
    # erfc(x) = 1 - erf(x)
    rho = 0.5 * erfc(arg) * (rho_A - rho_B) + rho_B
    return rho

def calc_f_N2(rho):
    return (rho_A * (rho_B - rho)) / (rho * (rho_B - rho_A))
    

def main():
    t = 1.0e-6
    fp = open("exact-profile.data", "w")
    fp.write("# x(m)    rho     f_N2     f_O2\n")
    x = xmin
    while( x < xmax ):
        rho = calc_rho( x, t )
        f_N2 = calc_f_N2( rho )
        f_O2 = 1.0 - f_N2
        fp.write("%.8f  %.8f  %.8f  %.8f  \n" % ( x, rho, f_N2, f_O2 ) )
        x += dx


if __name__ == '__main__': main()

