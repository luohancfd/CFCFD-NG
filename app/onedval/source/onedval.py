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
from prop_avg import stream_thrust_avg

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

def main():
    """
    Top-level function for the onedval program.
    """
    # 0. Gather command-line info
    if len(sys.argv) != 3:
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
    

    # 1a. Look at species and setup gas model
    if not 'species' in cfg:
        print "No 'species' list was found, so defaulting to ['air']."
        cfg['species'] = ['air']
    select_gas_model("thermally perfect gas", cfg['species'])
    gmodel = get_gas_model_ptr()

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

    # 2. Read data from slice
    cells = create_cells_from_slice(slice_file, cfg['variable_map'], cfg['grid_scale'])
    print "Total number of cells created from slice: ", len(cells)
    # 2a. apply filtering if required
    if 'filter_function' in cfg:
        print "Using filter function to remove unwanted or unimportant cells."
        cells = filter(cfg['filter_function'], cells)
        print "Number of cells after filtering: ", len(cells)

    
    # 3. Do some work.
    # 3a. Compute requested one_d_properties
    # -- TEST ---
    stream_thrust_avg(cells, [], cfg['variable_map'], cfg['species'], gmodel)
    

        
        
        
        

if __name__ == '__main__':
    main()
    

    
