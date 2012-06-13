# Author: Rowan J. Gollan
# Date: 02-Apr-2012
#
# A script to plot the Cp of diatomic oxygen
# over the temperature range of 200.0 -- 20,000.0 K.
#

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from gaspy import *
import numpy as np
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

# Initialise gas model and gas data
gmodel = create_gas_file('thermally perfect gas', ['O2'],
                'thermally-perfect-O2.lua')
gd = Gas_data(gmodel)

gd.massf[0] = 1.0
gd.p = 1.0e5 # Some sensible value to keep the module happy
             # does not affect the calculation of Cp

# Select the temperature values for Cp evaluation
start = 200.0
stop = 20000.0
step = 100.0
T_list = np.arange(start, stop + 0.5*step, step)

# Begin looping to evaulate Cp
Cp_list = np.zeros_like(T_list)
for i, T in enumerate(T_list):
    gd.T[0] = T
    Cp_list[i] = gmodel.Cp(gd)

# Plot result and exit
plt.axis([200.0, 20000.0, 0.0, 1500.0])
plt.xlabel('Temperature, K')
plt.ylabel(r'$C_p$, J/(kg.K)')
plt.plot(T_list, Cp_list)
plt.savefig('../figs/Cp-O2.pdf')

