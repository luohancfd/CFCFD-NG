#!/usr/bin/env python
# nenzfr_check_forebody.py
#
# A short script to compare the calculation of the post-forebody
# shock flow properties using a thermally-perfect gas model and
# and ideal gas model.
#
#
# Luke Doherty
# July-2012 ?

import shlex, string
import sys, os
import optparse, copy, math
from numpy import array, sqrt
from nenzfr_utils import run_command, quote, read_case_summary, \
     read_nenzfr_outfile, read_estcj_outfile, \
     read_gmodelFile_from_config
from libprep3 import create_gas_model, Gas_data, set_massf
import cfpylib.gasdyn.gas_flow as gf
from cfpylib.nm.zero_solvers import secant
from cfpylib.gasdyn.ideal_gas import Gas as Gas_ideal
from nenzfr_forebody_sensitivity import Gas as Gas_thermal


E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#
# This script must be called from the following directory:
#     work/space_planes_2012/sensitivity_air/
# 


# Load some nenzfr data for us to use
data, Var = read_nenzfr_outfile('./case000/nozzle-exit.stats')

# Read gas model file
gmodelFile = read_gmodelFile_from_config('./case000/nozzle')
if not os.path.exists(gmodelFile):
    run_command('cp ./case000/'+gmodelFile+' ./')

# Define forebody angle
theta_rad = 6.0*math.pi/180.0

# Extract a species dictionary from the input data
speciesKeys = [k for k in data.keys() if k.startswith("mass")]
speciesMassFrac = [data[k] for k in speciesKeys]
speciesDict = dict([(k.split('-')[1],v) for k,v in zip(speciesKeys,speciesMassFrac)])

# Create a Thermally Perfect Gas object
thermal = Gas_thermal(name='air5species',speciesDict=speciesDict,gasModelFile=gmodelFile)
thermal.set_pT(p=data['p'],T=data['T[0]'])

# Create an Ideal Gas object
ideal = Gas_ideal(Mmass=8314.0/thermal.R, gamma=thermal.gam, name='air5species')
ideal.set_pT(p=data['p'],T=data['T[0]'])

# Define velocity
vel = sqrt(data['vel.x']*data['vel.x'] + data['vel.y']*data['vel.y'])

# Now calculate post-forebody shock conditions for thermally perfect gas
thermal_result = {}
thermal_result['beta'] = gf.beta_oblique(thermal, vel, theta_rad)

thermal_result['theta_out'], thermal_result['fbV'], thermalFBstate = \
                              gf.theta_oblique(thermal, vel, thermal_result['beta'])

thermalPitot = gf.pitot_condition(thermal, vel)
thermalTotal = gf.total_condition(thermal, vel)

print "THERMAL RESULT:"
print thermal_result
thermalFBstate.write_state(sys.stdout)
print "pitot:"
thermalPitot.write_state(sys.stdout)
print "total:"
thermalTotal.write_state(sys.stdout)


# Calculate post-forebody shock conditions for ideal gas
ideal_result = {}
ideal_result['beta'] = gf.beta_oblique(ideal, vel, theta_rad)

ideal_result['theta_out'], ideal_result['fbV'], idealFBstate = \
                             gf.theta_oblique(ideal, vel, ideal_result['beta'])

idealPitot = gf.pitot_condition(ideal, vel)
idealTotal = gf.total_condition(ideal, vel)

print 
print "IDEAL RESULT:"
print ideal_result
idealFBstate.write_state(sys.stdout)
print "pitot:"
idealPitot.write_state(sys.stdout)
print "total:"
idealTotal.write_state(sys.stdout)


