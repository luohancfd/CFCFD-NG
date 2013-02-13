#!/usr/bin/env python
"""
Python program to compute one-dimensionalized quantities from 2D data.

This program is supplied with two inputs by the user:
1) a previously prepared Tecplot slice of data in ASCII block format; and
2) a config file to coordinate the calculation.
The program will then compute one-dimensionalized quantities from the
the 2D surface based on the methods selected by the user.

.. Author: Rowan J. Gollan

.. Versions:
   12-Feb-2013  Initial coding.
"""

import sys

from e3prep import select_gas_model
from libprep3 import get_gas_model_ptr
from cell import create_cells_from_slice
from prop_avg import *

default_var_map = {'x':'x', 'y':'y', 'z':'z', 'u':'u', 'v':'v', 'w':'w',
                   'rho':'rho', 'p':'p', 'T':'T', 'M':'M', 'h0':'h0'}

def print_usage():
    print ""
    print "Usage: onedval INPUT CONFIG [OUTPUT]"
    print ""
    print "where:"
    print "INPUT  -- name of Tecplot file with slice data"
    print "CONFIG -- name of config file to control calculation"
    print "OUTPUT -- optional name of file for output data"
    print "          if not supplied default output is 'onedval-props.txt'"
    print ""

pretty_var_names = {'rho':'density (kg/m^3)',
                    'p':'pressure (Pa)',
                    'T':'temperature (K)',
                    'M':'Mach number',
                    'u':'U velocity (m/s)',
                    'v':'V velocity (m/s)',
                    'w':'W velocity (m/s)'}

def pretty_print_props(f, props, species):
    for k,v in props.iteritems():
        if k in pretty_var_names:
            f.write("%s\n" % pretty_var_names[k])
            f.write("%s = %.6e\n" % (k, v))
        elif k in species:
            f.write("mass fraction of %s : %.6e\n" % (k, v))
        else:
            f.write("%s : %.6e\n" % (k, v))
        f.write("\n")
    return

def main():
    """
    Top-level function for the onedval program.
    """
    print "onedval: A program to compute integrated and one-dimensionalised quantities."
    print "onedval: Beginning."
    # 0. Gather command-line info
    if len(sys.argv) < 3:
        print "At least two arguments are required."
        print_usage()
        sys.exit(1)

    slice_file = sys.argv[1]
    config_file = sys.argv[2]
    output_file = 'onedval-props.txt'
    if len(sys.argv) == 4:
        output_file = sys.argv[3]
    
    # 1. Gather info from config file
    # Set some defaults.
    # If set to 'None', we expect to find something from the user.
    cfg = {}
    try:
        execfile(config_file, globals(), cfg)
    except IOError:
        print "There was a problem reading the config file: ", config_file
        print "Check that is conforms to Python syntax."
        print "Bailing out!"
        sys.exit(1)
    
    print "onedval: Setting up gas model"
    # 1a. Look at species and setup gas model
    if not 'species' in cfg:
        print "No 'species' list was found, so defaulting to ['air']."
        cfg['species'] = ['air']
    select_gas_model("thermally perfect gas", cfg['species'])
    gmodel = get_gas_model_ptr()
    nsp = gmodel.get_number_of_species()

    print "onedval: Checking over user inputs"
    # 1b. Check for variable map
    if not 'variable_map' in cfg:
        print "No 'variable_map' was set so the following default is used:"
        print default_var_map
        cfg['variable_map'] = default_var_map
    
    # 1c. Look for one_d_averages methods
    if not 'one_d_averages' in cfg:
        print "No list of 'one_d_averages' was set."
        print "The default method of 'flux-conserved' will be used."
        cfg['one_d_averages'] = ['flux-conserved']
    
    # 1d. Look for grid_scale
    if not 'grid_scale' in cfg:
        print "No 'grid_scale' was set."
        print "The default value of 1.0 will be used."
        cfg['grid_scale'] = 1.0

    # 1e. Look for one_d_outpus
    if not 'one_d_outputs' in cfg:
        cfg['one_d_outputs'] = []

    # 1f. Look for integrated_outputs
    if not 'integrated_outputs' in cfg:
        cfg['integrated_outputs'] = []

    print "onedval: Reading in data from slice and creating cells"
    # 2. Read data from slice
    cells = create_cells_from_slice(slice_file, cfg['variable_map'], cfg['grid_scale'])
    print "Total number of cells created from slice: ", len(cells)
    # 2a. apply filtering if required
    if 'filter_function' in cfg:
        print "Using filter function to remove unwanted or unimportant cells."
        cells = filter(cfg['filter_function'], cells)
        print "Number of cells after filtering: ", len(cells)

    
    print "onedval: Doing the requested calculations"
    # 3. Do some work.
    f = open(output_file, 'w')
    f.write("------------------- onedval output ---------------------\n\n")
    f.write("number of cells in averaging:\n")
    f.write("ncells = %d\n" % len(cells))
    f.write("cumulative area of cells (m^2):\n")
    f.write("area = %.6e\n" % area(cells))

    # 3a. Compute requested integrated quantities (if required)
    if len(cfg['integrated_outputs']) > 0:
        print "onedval: Writing out integrated quantities"
        f.write("\n---------------------\n")
        f.write("Integrated quantities\n")
        f.write("---------------------\n")
        f.write("\n")

    fluxes = compute_fluxes(cells, cfg['variable_map'], cfg['species'], gmodel)

    for flux in cfg['integrated_outputs']:
        if flux == 'mass flux':
            f.write("mass flux (kg/s)\n")
            f.write("m_dot = %.6e\n" % fluxes['mass'])
        elif flux == 'momentum flux':
            f.write("momentum flux (kg.m/s^2)\n")
            f.write("mom_dot = %.6e\n" % fluxes['mom'])
        elif flux == 'energy flux':
            f.write("energy flux (W)\n")
            f.write("e_dot = %.6e\n" % fluxes['energy'])
        elif flux == 'species mass flux':
            if nsp == 1:
                continue
            for sp in cfg['species']:
                isp = gmodel.get_isp_from_species_name(sp)
                f.write("mass flux of %s (kg/s)\n" % sp)
                f.write("m%s_dot = %.6e\n" % (sp, fluxes['species'][isp]))
        else:
            print "Requested integrated quantity: ", f
            print "is not part of the list of available integrated quantities."
            print "Bailing out!"
            sys.exit(1)

    # 3b. Compute requested one_d_properties
    if len(cfg['one_d_averages']) > 0:
        print "onedval: Writing out one-dimensionalised quantities"
        f.write("\n------------------------------\n")
        f.write("One-dimensionalised quantities\n")
        f.write("------------------------------\n")
        f.write("\n")


    for avg in cfg['one_d_averages']:
        if avg == 'area-weighted':
            f.write("=== area-weighted average\n\n")
            phis = area_weighted_avg(cells, cfg['one_d_outputs'], cfg['variable_map'])
            pretty_print_props(f, phis, cfg['species'])
            f.write("\n")
        elif avg == 'mass-flux-weighted':
            f.write("=== mass flux weighted average\n\n")
            phis = mass_flux_weighted_avg(cells, cfg['one_d_outputs'], cfg['variable_map'])
            pretty_print_props(f, phis, cfg['species'])
            f.write("\n")
        elif avg == 'flux-conserved':
            f.write("=== flux-conserved average\n\n")
            phis = stream_thrust_avg(cells, cfg['one_d_outputs'], cfg['variable_map'], cfg['species'], gmodel)
            pretty_print_props(f, phis, cfg['species'])
        else:
            print "Requested one-D averaging method: ", avg
            print "is not known or not implemented."
            print "Bailing out!"
            sys.exit(1)

    f.close()
    print "onedval: Done."

if __name__ == '__main__':
    main()
    

    
