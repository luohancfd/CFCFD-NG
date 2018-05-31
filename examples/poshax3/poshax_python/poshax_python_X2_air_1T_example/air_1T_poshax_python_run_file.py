#! /usr/bin/env python
"""
air_1T_poshax_python_run_file:

This is a one temperature air example for the poshax_python program
using some proposed X2 conditions.

A two temperature version of this example is also available.

Chris James (c.james4@uq.edu.au) - 27/04/18

"""

import sys, os, math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from poshax_python import poshax_python

# these are freestream / nozzle exit conditions for a proposed X2 expansion tube condition
# where we were trying to achieve a longer post shock relaxation condition
title = "Steve Apirana's random condition" # title for the plots
p_inf = 147.72 # Pa
T_inf = 2242.5 # K
u_inf = 9690.6 # m/s

reactants = {'N2': 0.79, 'O2': 0.21}  # reactant dict for cea calc
inputUnits = 'moles'  # for cea
outputUnits = 'massf' # for poshax output...

species_list = ['N2', 'N2_plus', 'NO', 'NO_plus', 'O2', 'O2_plus', 'N', 'N_plus', 'O', 'O_plus', 'e_minus']

no_of_temperatures = 1

# name of the .inp gas file which we will make...
# and the .lua gas file which the gasfile program will make
# don't add the .lua here as we have two different files to deal with here
gas_model_file = 'air-11sp-1T'

# these are both lua files so add the .lua
# reaction file is always needed
reaction_file = 'gupta_etal_air_reactions.lua'
# energy exchange file is only needed for multiple temperatures
energy_exchange_file = None

# cfg file is the file whcih we'll send to poshax to run the simulation
cfg_file = 'cfg_file.cfg'

# output_file is the result table with distance from the shock
output_file = 'output_file.data'

# simulation values for poshax
source_term_coupling = 'loose'
dx = 1.0e-16
adaptive_dx = 'true'
final_x = 5.0e-3

plot = True
save_figures = True
log_x_axis = False
plot_x_limits= [0.0,final_x*1.0e3] # mm

poshax_python(p_inf=p_inf, T_inf=T_inf, u_inf=u_inf,
              reactants=reactants, inputUnits=inputUnits, outputUnits=outputUnits,
              species_list=species_list, no_of_temperatures=no_of_temperatures,
              reaction_file=reaction_file, energy_exchange_file=energy_exchange_file,
              gas_model_file=gas_model_file, cfg_file=cfg_file, output_file = output_file,
              source_term_coupling=source_term_coupling, dx=dx, adaptive_dx=adaptive_dx,
              final_x=final_x, plot=plot, save_figures=save_figures, title=title,
              log_x_axis=log_x_axis, plot_x_limits=plot_x_limits)
