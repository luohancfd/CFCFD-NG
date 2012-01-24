# test_noble_abel_gas.py
# PJ 08-feb-2011

import sys
import os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from gaspy import *

def percent_diff(x, y):
    return 100.0 * abs(x-y) / max(abs(x),abs(y))

print "Begin test of Noble-Abel gas model."
p = 500.0e6 # Pa
T = 2000.0  # degree K
b = 0.001
print "p=", p, "T=", T, "b=", b

print "\nNoble-Abel gas:"
gmodel = create_gas_model('ja2-propellant-gas.lua')
gas = Gas_data(gmodel)
gas.p = p; gas.T[0] = T; gas.massf[0] = 1.0
gmodel.eval_thermo_state_pT(gas)
gmodel.eval_transport_coefficients(gas)
gas.print_values()
R = gmodel.R(gas); Cv = gmodel.Cv(gas); Cp = gmodel.Cp(gas); gamma = gmodel.gamma(gas)
print "R=", R, "Cv=", Cv, "Cp=", Cp, "gamma=", gamma
rho = p / (R*T + b*p)
print "expected density=", rho, "percent-difference=", percent_diff(rho, gas.rho)
e = Cv * T
print "expected internal energy=", e, "percent-difference=", percent_diff(e, gas.e[0])

print "\nIdeal gas:"
create_gas_file(model='ideal gas', species=['air'], fname='ideal-propellant.lua')
change_ideal_gas_attribute('air', 'M', 0.02489, fname="ideal-propellant.lua")
change_ideal_gas_attribute('air', 'gamma', 1.225, fname="ideal-propellant.lua")
gmodel_ideal = create_gas_model('ideal-propellant.lua')
gas2 = Gas_data(gmodel)
gas2.p = p; gas2.T[0] = T; gas2.massf[0] = 1.0
gmodel_ideal.eval_thermo_state_pT(gas2)
gmodel_ideal.eval_transport_coefficients(gas2)
gas2.print_values()

print "\nDifferences with respect to Noble-Abel gas:"
print "density=", percent_diff(gas2.rho, gas.rho)
print "sound speed=", percent_diff(gas2.a, gas.a)
print "internal energy=", percent_diff(gas2.e[0], gas.e[0])

print "Done."
