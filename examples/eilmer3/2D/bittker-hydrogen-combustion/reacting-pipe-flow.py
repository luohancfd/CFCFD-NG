# reacting-pipe-flow.py
# PJ 21-Jun-2011
#
# Combustion in a supersonic stream of constant cross-section.
# This is set up to approximate the Bittker-Scullin case 3
# as used by Fabs in the hydrogen-combustion test case. 

#-------------------------------------------------------------
# Things that we'll make use of shortly.

import sys, os, math
import numpy as np
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from libprep3 import *

def sample_header():
    return "# x(m) rho(kg/m**3) p(Pa) T(degK) e(J/kg) v(m/s) "\
        "massf_OH massf_H2O dt_suggest(s)"
def sample_data(x, v, gas, dt_suggest):
    return "%g %g %g %g %g %g %g %g %g " % (x, gas.rho, gas.p, gas.T[0], gas.e[0],
                                           v, gas.massf[7], gas.massf[5], dt_suggest)

def eos_derivatives(gas0, gmodel, tol=0.0001):
    """
    Finite difference evaluation, assuming that gas0 is valid state already.
    """
    gas1 = Gas_data(gas0)
    p0 = gas0.p
    rho0 = gas0.rho
    e0 = gas0.e[0]
    #
    drho = rho0 * tol
    gas1.rho = rho0 + drho
    gmodel.eval_thermo_state_rhoe(gas1)
    dpdrho = (gas1.p - p0)/drho
    #
    gas1.rho = rho0
    de = e0 * tol
    gas1.e[0] = e0 + de
    gmodel.eval_thermo_state_rhoe(gas1)
    dpde = (gas1.p - p0)/de
    #
    return dpdrho, dpde

#------------------------------------------------------------------
# Start the main script...
debug = False
do_gas_dynamic_accommodation = True

species = ['O', 'O2', 'N2', 'H', 'H2', 'H2O', 'HO2', 'OH', 'H2O2']
create_gas_file('thermally perfect gas', species, 'h2-air.lua')
gmodel = create_gas_model('h2-air.lua')
nsp = gmodel.get_number_of_species()

print "# Reacting pipe flow -- Bittker-Scullin test case 3." 
print sample_header()
print "# Gas properties at the start of the pipe."
molef = {'O2':0.1480, 'N2':0.5562, 'H2':0.2958}
massf = gmodel.to_massf(molef)
gas0 = Gas_data(gmodel)
gas0.p = 96.87e3 # Pa
x = 0.0 # m  (inlet of pipe)
v = 4551.73 # m/s
gas0.T[0] = 1559.0 # degree K
for i in range(nsp): gas0.massf[i] = massf[i]
gmodel.eval_thermo_state_pT(gas0)
dt_suggest = 1.0e-8  # suggested starting time-step for chemistry updater
print sample_data(x, v, gas0, dt_suggest)

print "# Start reactions..."
rupdate = create_Reaction_update('Bittker_Scullin.lua', gmodel)
t = 0 # time is in seconds
t_final = 22.0e-6
t_inc = 0.1e-6
nsteps = int(t_final / t_inc)
for j in range(nsteps):
    # At the start of the step...
    rho = gas0.rho
    T = gas0.T[0]
    p = gas0.p
    e = gas0.e[0]
    #
    # Do the chemical increment.
    gas1 = Gas_data(gas0) # make the new one as a clone
    dt_suggest = rupdate.update_state_py(gas1, t_inc, dt_suggest, gmodel)
    gmodel.eval_thermo_state_rhoe(gas1)
    # gmodel.eval_transport_coefficients(gas1)
    #
    de_chem = gas1.e[0] - e
    dp_chem = gas1.p - p
    if debug: print "# de_chem=", de_chem, "dp_chem=", dp_chem
    #
    if do_gas_dynamic_accommodation:
        # Do the gas-dynamic accommodation after the chemical change.
        Etot = e + 0.5*v*v
        dfdr, dfde = eos_derivatives(gas1, gmodel)
        if debug: print "# dfdr=", dfdr, "dfde=", dfde
        A = np.array([[v,      rho,        0.0, 0.0  ],
                      [0.0,    rho*v,      1.0, 0.0  ],
                      [v*Etot, rho*Etot+p, 0.0, rho*v],
                      [dfdr,   0.0,       -1.0, dfde ]]);
        b = np.array([0.0, -dp_chem, -rho*v*de_chem, 0.0])
        dq = np.linalg.solve(A,b)
        drho, dv, dp_gda, de_gda = dq
        if debug: 
            print "# A*dq-b=", np.dot(A,dq)-b
            print "# drho=", drho, "dv=", dv, "dp_gda=", dp_gda, "de_gda=", de_gda
        # Add the accommodation increments.
        gas1.rho = gas0.rho + drho
        v1 = v + dv
        p1_check = gas1.p + dp_gda
        gas1.e[0] += de_gda
        gmodel.eval_thermo_state_rhoe(gas1)
        if debug: print "# At new point: gas1.p=", gas1.p, "p1_check=", p1_check, \
            "rel_error=", abs(gas1.p-p1_check)/p1_check
    else:
        # Don't do the gas-dynamic accommodation.
        v1 = v
    # Have now finished the chemical and gas-dynamic update.
    t += t_inc
    x += 0.5*(v + v1) * t_inc
    print sample_data(x, v1, gas1, dt_suggest)
    # House-keeping for the next step.
    v = v1
    gas0.copy_values_from(gas1)

print "# Done."
