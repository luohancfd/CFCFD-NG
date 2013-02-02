#!/usr/bin/env python
# nenzfr_eq_cone_pitot_pressures.py
#
# This script loads in a nominated nenzfr slice data
# file and calculates the expected profile for the 
# cone static pressure and pitot pressure assuming an
# equilibrium gas. The profiles are saved to a .data 
# file.
# 
# Luke Doherty
# School of Mechanical and Mining Engineering
# The University of Queensland

VERSION_STRING = "10-July-2012"

import sys, os, optparse, string, numpy
from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.gas_flow import *

E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

def load_nenzfr_exit_data(dataFile):
    """
    Load in a nenzfr slice data file
    """
    fp = open(dataFile, 'r')
    # Keep a list of variables in order of appearance.
    varLine = fp.readline().strip()
    items = varLine.split()
    if items[0] == '#': del items[0]
    if items[0] == 'Variables:': del items[0]
    variable_list = [item.split(':')[1] for item in items]
    # print "variable_list=", variable_list
    # Store the data in lists against these names.
    data = {}
    for var in variable_list:
        data[var] = []
    for line in fp.readlines():
        items = line.strip().split()
        if items[0] == '#': continue
        assert len(items) == len(variable_list)
        for i in range(len(items)):
            data[variable_list[i]].append(float(items[i]))
    fp.close()
    
    return data, variable_list

def calculate_pressure(data, variable_list, theta):
    """
    Using the input data, determine the gas properties
    and then calculate the expected pressure on the 
    surface of a cone of half-angle theta assuming an
    equilibrium gas.
    """ 

    theta_rad = theta*numpy.pi/180.0;
    
    # Find the species present in data set
    species_onlyList = []
    for name in variable_list:
        if name.startswith('massf'):
            speciesName = name.split('-')[-1]
            species_onlyList.append(speciesName)
    # Initialise some variables
    eq_pitot_p = []; cone_beta = []; cone_theta = []; 
    cone_V = []; cone_p = []; cone_T = [];

    for j in range(len(data['pos.y'])):
        # Create Gas object for current location...
        react = {}
        for spc in species_onlyList:
            react[spc] = data['massf['+\
                    str(species_onlyList.index(spc))+']-'+spc][j]
        #print "reactants=",react
        gasState = Gas(reactants=react, inputUnits='massf',\
                       onlyList=species_onlyList,\
                       outputUnits='massf')
        # Set pressure and temperature...
        gasState.set_pT(p=data['p'][j], T=data['T[0]'][j])
        
        # Calculate the pitot pressure
        pitotState = pitot_condition(gasState, data['vel.x'][j])
        eq_pitot_p.append( pitotState.p )
        
        
        #print eq_pitot_p
        
        # Calculate the cone static pressure
        #
        # Only if the Mach number is greater than 1 do we calculate
        # the surface pressure on the cone using the Taylor-Maccoll
        # equations. Otherwise we just set it equal to the freestream
        # pressure (we aren't interested in these low supersonic
        # regions).
        if data['M_local'][j] > 2.1: #1.05:
            # Calculate the conical shock angle
            cone_beta.append( beta_cone(gasState, data['vel.x'][j], theta_rad) )
            #print cone_beta
            # Calculate the gas properties on the cone surface
            theta_c, V_c, coneState = theta_cone(gasState, data['vel.x'][j], cone_beta[j])
        
            cone_theta.append( theta_c )
            cone_V.append( V_c )
            cone_p.append( coneState.p )
            cone_T.append( coneState.T )
            
            #assert abs(theta_c-theta_rad) < 1.0e-4
            print "theta_error=",abs(theta_c-theta_rad)
        else:
            cone_p.append( data['p'][j] )
            cone_T.append( data['T[0]'][j] )
            cone_V.append( 0.0 )
            cone_theta.append( 0.0 )
            cone_beta.append( 0.0 )
        
        print "j = {0:d} of {1:d} complete".format(j,len(data['pos.y']))
        
    variable_list.append('cone_p'); variable_list.append('cone_T')
    variable_list.append('cone_beta'); variable_list.append('cone_theta')
    variable_list.append('cone_V'); variable_list.append('eq_pitot_p')
    
    data['cone_p'] = cone_p; data['cone_T'] = cone_T;
    data['cone_beta'] = cone_beta; data['cone_theta'] = cone_theta; 
    data['cone_V'] = cone_V;
    data['eq_pitot_p'] = eq_pitot_p;
    
    return data, variable_list
    

def write_output(data, variable_list, outFileName):
    """
    """
    fout = open(outFileName,'w')
    fout.write('# File produced by: nenzfr_eq_cone_pitot_pressure.py\n')
    fout.write('# Variables: 1:pos.x 2:pos.y 3:cone_p 4:cone_T 5:cone_beta 6:cone_theta 7:eq_pitot_p \n')
    for j in range(len(data['pos.x'])):
        fout.write('{0:7.6e} {1:7.6e} '.format(data['pos.x'][j], data['pos.y'][j]))
        fout.write('{0:7.6e} {1:7.6e} '.format(data['cone_p'][j], data['cone_T'][j]))
        fout.write('{0:7.6e} {1:7.6e} '.format(data['cone_beta'][j],data['cone_theta'][j]))
        fout.write('{0:7.6e} '.format(data['eq_pitot_p'][j]))
        fout.write('\n')
    fout.close()

def main():
    """
    """
    op = optparse.OptionParser(version=VERSION_STRING)
    
    op.add_option('--run-defaults', dest='runDefaults', action='store_true',
                  default=True, help="calculate properties on surface "
                  "of a cone using all default paramters.")
    op.add_option('--datafile', dest='exitDataFile',
                  default='nozzle-exit.data',
                  help="file that holds the nozzle-exit data "
                       "[default: %default]")
    op.add_option('--outfile', dest='coneDataFile',
                  default='nozzle-exit-cone-pitot.data',
                  help="file to write the cone surface and pitot data to "
                       "[default: %default]")
    op.add_option('--theta', dest='theta',type=float, default=6.0,
                  help="Cone half angle in degrees [default: %default]")
    opt, args = op.parse_args()
    
    # Load data
    data, variable_list = load_nenzfr_exit_data(opt.exitDataFile)
    
    # Calculate pressure on surface of a cone
    data, variable_list = \
          calculate_pressure(data, variable_list, opt.theta)
    
    # Write the data out to a file suitable for use with gnuplot
    write_output(data, variable_list, opt.coneDataFile)


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "nenzfr_eq_cone_pitot_pressure: Calculate the expected cone static and pitot pressure profiles across the nozzle exit assuming an equilibrium gas."
        print "    Version:", VERSION_STRING
        print "    To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
