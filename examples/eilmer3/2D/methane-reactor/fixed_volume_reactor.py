#!/usr/bin/env python
"""
Simple fixed-volume reactor

This script uses the thermochemical module directly and
replaces Brendan's more complex reactor code.

Peter J. 24-Apr-2013
"""

import sys
import os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from gaspy import *
from numpy import arange

gmodel = create_gas_model('thermally-perfect-grimech30.lua')
Q = Gas_data(gmodel)
set_molef(Q, gmodel, {'CH4':1.0, 'O2':2.0, 'N2':7.52})
Q.T[0] = 2000.0 # degrees K
Q.p = 101.325e3 # Pascals
gmodel.eval_thermo_state_pT(Q)
nsp = gmodel.get_number_of_species()
nmodes = gmodel.get_number_of_modes()
if 0:
    print "Initial gas state: p=", Q.p, "rho=", Q.rho
    print "    massf=", [Q.massf[i] for i in range(nsp)]
    print "    T=", [Q.T[i] for i in range(nmodes)]
    print "    e=", [Q.e[i] for i in range(nmodes)]

def write_header():
    """
    Write a header in the same format as Brendan's program.
    """
    print "# t, T_0, p_0, T, p, rho",
    for i in range(nsp): print (", '%s'" % gmodel.species_name(i)),
    print
    return
def write_data(tme):
    """
    Write gas data in the same format as Brendan's program.
    """
    print ("%.5e %.5e %.5e %.5e %.5e %.12e" % (tme, Q.T[0], Q.p, Q.T[0], Q.p, Q.rho)),
    molef = convert_massf2molef(Q.massf, gmodel.M())
    for i in range(nsp):
        print (" %.12e" % molef[i]),
    print
    return

r = create_Reaction_update('grimech30.lua', gmodel)
t_end = 4.0e-4 # seconds
dt = t_end/2000
dt_chem = -1.0 # Let Rowan's module choose
all_times = arange(0.0, t_end, dt)
write_header()
write_data(all_times[0])


print "# Now, integrate in time."
for t in all_times[1:]:
    dt_chem = r.update_state_py(Q, dt, dt_chem, gmodel)
    gmodel.eval_thermo_state_rhoe(Q)
    write_data(t)

print "# Done."
 
