#!/usr/bin/env python
# nenzfr_cone_static_pressure.py
#
# This script loads in a nominated nenzfr slice data
# file and calculates the expected profile for the 
# cone static pressure thereby allowing direct 
# comparison to experimental survey data.
#
# It assumes and IDEAL GAS
# 
# Luke Doherty
# School of Mechanical and Mining Engineering
# The University of Queensland

VERSION_STRING = "02-July-2012"

import sys, os, optparse, string, numpy
from cfpylib.gasdyn.ideal_gas_flow import taylor_maccoll_odes,\
    theta_cone, beta_cone 
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

def calculate_cone_pressure(data, variable_list, theta):
    """
    Using the input data, determine the gas properties
    and then calculate the expected pressure on the 
    surface of a cone of half-angle theta assuming an
    ideal gas.
    """ 

    theta_rad = theta*numpy.pi/180.0;
    
    R = []; g = []; beta = []; theta_c = []; V_r = [];
    p_c = []; T_c = [];

    for j in range(len(data['pos.y'])):
        # Calculate Gas Constant
        R.append( data['p'][j]/(data['rho'][j]*data['T[0]'][j]) )
        # Calculate Ratio of Specific Heats
        g.append( data['a'][j]*data['a'][j]/(R[j]*data['T[0]'][j]) )
        
        # Only if the Mach number is greater than 1 do we calculate
        # the surface pressure on the cone using the Taylor-Maccoll
        # equations. Otherwise we just set it equal to the freestream
        # pressure (we aren't interested in these low supersonic
        # regions).
        if data['M_local'][j] > 1.05:
            # Calculate conical shock angle
            beta.append( beta_cone(data['vel.x'][j],data['p'][j],\
                                   data['T[0]'][j],theta_rad,\
                                   R[j], g[j]) )
            # Calculate the gas properties on the cone surface
            surface = theta_cone(data['vel.x'][j],data['p'][j],\
                                 data['T[0]'][j],beta[j],\
                                 R[j], g[j])
            theta_c.append( surface[0] )
            V_r.append( surface[1] )
            p_c.append( surface[2] )
            T_c.append( surface[3] )
            
            assert abs(theta_c[j]-theta_rad) < 1.0e-4
        else:
            p_c.append( data['p'][j] )
            T_c.append( data['T[0]'][j] )
            V_r.append( 0.0 )
            theta_c.append( 0.0 )
            beta.append( 0.0 )
            
    variable_list.append('p_c'); variable_list.append('T_c')
    variable_list.append('beta_c'); variable_list.append('theta_c')
    variable_list.append('V_r'); variable_list.append('R')
    variable_list.append('gamma')
    
    data['p_c'] = p_c; data['T_c'] = T_c; data['R'] = R; data['gamma'] = g;
    data['beta_c'] = beta; data['theta_c'] = theta_c; data['V_r'] = V_r;
    
    return data, variable_list
    

def write_output(data, variable_list, outFileName):
    """
    """
    fout = open(outFileName,'w')
    fout.write('# Variables: 1:pos.x 2:pos.y 3:p_c 4:T_c 5:beta_c 6:theta_c 7:R 8:gamma \n')
    for j in range(len(data['pos.x'])):
        fout.write('{0:7.6e} {1:7.6e} '.format(data['pos.x'][j], data['pos.y'][j]))
        fout.write('{0:7.6e} {1:7.6e} '.format(data['p_c'][j], data['T_c'][j]))
        fout.write('{0:7.6e} {1:7.6e} '.format(data['beta_c'][j],data['theta_c'][j]))
        fout.write('{0:7.6e} {1:7.6e} '.format(data['R'][j], data['gamma'][j]))

        #data['pos.x'][j], data['pos.y'][j], \
        #data['p_c'][j], data['T_c'][j], data['beta_c'][j],\
        #data['theta_c'][j], data['R'][j], data['gamma'][j])
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
                  default='nozzle-exit-cone-surface.data',
                  help="file to write the cone surface data to "
                       "[default: %default]")
    op.add_option('--theta', dest='theta',type=float, default=6.0,
                  help="Cone half angle in degrees [default: %default]")
    opt, args = op.parse_args()
    
    # Load data
    data, variable_list = load_nenzfr_exit_data(opt.exitDataFile)

    # Calculate pressure on surface of a cone
    data, variable_list = \
          calculate_cone_pressure(data, variable_list, opt.theta)
    
    # Write the data out to a file suitable for use with gnuplot
    write_output(data, variable_list, opt.coneDataFile)


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "nenzfr_cone_static_pressure: Calculate the expected cone static pressure profile across the nozzle exit."
        print "    Version:", VERSION_STRING
        print "    To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
