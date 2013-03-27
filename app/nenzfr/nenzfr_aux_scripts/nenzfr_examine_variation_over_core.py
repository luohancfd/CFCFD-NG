#!/usr/bin/env python
# nenzfr_examine_variation_over_core.py
#
# This script loads in data from the cases used to build 
# a response surface and then writes out a file summarising
# the standard deviation for the variation over the core
# for each flow property and each simulation case. 
#
# Luke Doherty
# School of Mechanical and Mining Engineering
# The University of Queensland

VERSION_STRING = "27-March-2013"

import string
import sys, os
import optparse
import numpy as np
from nenzfr_utils import read_case_summary, read_nenzfr_outfile
from nenzfr_stats import get_slice_data
#from nenzfr_sensitivity import add_extra_variables
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

def add_extra_variables(data, var):
    """
    We add some additional variables that of are interest to
    the input dictionary and list.

    data: dictionary of data
    var: list of variables
    """
    U_squared = data['vel.x']['mean']*data['vel.x']['mean'] + \
                data['vel.y']['mean']*data['vel.y']['mean']
    U = U_squared**0.5
    # Dynamic Pressure
    data['q'] = {'mean':0.5 * data['rho']['mean'] * U_squared,
                 'min':0.0, 'max':0.0, 'std':0.0}
    var.append('q')
    # Mass flow rate per unit area
    data['m_dot'] = {'mean':data['rho']['mean'] * U,
                     'min':0.0, 'max':0.0, 'std':0.0}
    var.append('m_dot')
    # Unit Reynolds number
    data['Re_u'] = {'mean':data['m_dot']['mean'] / data['mu']['mean'],
                    'min':0.0, 'max':0.0, 'std':0.0}
    var.append('Re_u')
    # Pressure coefficient
    data['p/q'] = {'mean':data['p']['mean'] / data['q']['mean'],
                   'min':0.0, 'max':0.0, 'std':0.0}
    var.append('p/q')

    return data, var

def add_extra_variables_to_slice(data, var):
    """
    Copied from 'nenzfr_sensitivity.py'
    """
    U_squared = np.array(data['vel.x'])*np.array(data['vel.x']) +\
                 np.array(data['vel.y'])*np.array(data['vel.y'])
    U = U_squared**0.5
    # Dynamic Pressure
    data['q'] = 0.5 * np.array(data['rho']) * U_squared
    var.append('q')
    # Mass flow rate per unit area
    data['m_dot'] = np.array(data['rho']) * U
    var.append('m_dot')
    # Unit Reynolds number
    data['Re_u'] = np.array(data['m_dot']) / np.array(data['mu'])
    var.append('Re_u')
    # Pressure coefficient
    data['p/q'] = np.array(data['p']) / np.array(data['q'])
    var.append('p/q')

    return data, var


def main():
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('--exitSliceFile', dest='exitSliceFile',
                  default='nozzle-exit.data',
                  help=("file for holding the nozzle-exit data "
                       "[default: %default]"))

    op.add_option('--exitStatsFile', dest='exitStatsFileName',
                  default='nozzle-exit.stats',
                  help="file that holds the averaged nozzle-exit "
                       "data and is to be read in for each perturbation "
                       "case [default: %default]")

    op.add_option('--run-defaults', dest='runDefaults', action='store_true',
                  default=True, help="run code using all default inputs.")
    opt, args = op.parse_args()
    
    # Hard-code the core radius fraction because I'm too lazy to
    # write the necessary code to read it from the nenzfr_outfile
    coreRfraction = 0.317460
    
    # Read the "perturbation_cases.dat" file
    perturbedVariables, DictOfCases = read_case_summary()

    # Load in all the freestream data, including the min, max and std
    nozzleData = {}
    for case in DictOfCases.keys():
        nozzleData[case], exitVar = \
            read_nenzfr_outfile('./'+case+'/'+opt.exitStatsFileName,inclStats=1)
        # Add extra variables of interest (q, Re_u, m_dot, p/q)
        nozzleData[case], exitVar = \
            add_extra_variables(nozzleData[case], exitVar)
        
        # We have to do some stuffing around in order to add extra variables
        # of interest (q, Re_u, m_dot, p/q) and appropriately calculate the
        # max, min and standard deviation
        
        # Get the full set of slice data and add the extra variables to it
        var_list, sliceData = get_slice_data('./'+case+'/'+opt.exitSliceFile)
        sliceData, var_list = add_extra_variables_to_slice(sliceData, var_list)
        
        # Following is copied in part from "nenzfr_stats.py" 
        ys = sliceData['pos.y']
        y_edge = ys[-1] * coreRfraction
        
        for var in ['q', 'm_dot', 'Re_u', 'p/q']:
            # Identify low and high values.
            diff_minus = 0.0
            diff_plus = 0.0
            count = 0.0
            stddev = 0.0
            for j in range(len(ys)):
                if ys[j] > y_edge: break
                diff = sliceData[var][j] - nozzleData[case][var]['mean']
                diff_minus = min(diff, diff_minus)
                diff_plus = max(diff, diff_plus)
                count += 1
                stddev += diff**2
            # Calculate the sample standard deviation
            stddev = (stddev/(count-1))**0.5
        
            # Add data to the relevant dictionary
            nozzleData[case][var]['std'] = stddev
            nozzleData[case][var]['min'] = diff_minus
            nozzleData[case][var]['max'] = diff_plus
    
    # Calculate the standard-deviations as percentages of the 1D flow
    # value 
    std_as_percentage = {}
    for case in DictOfCases.keys():
        std_as_percentage[case] = {}
        for pr in exitVar:
            if nozzleData[case][pr]['mean'] != 0.0:
                std_as_percentage[case][pr] = nozzleData[case][pr]['std']/\
                                              nozzleData[case][pr]['mean']\
                                              * 100.0
            else:
                std_as_percentage[case][pr] = float('NaN')
    
    # Write an output file
    fp = open('summary_of_variation_across_core.dat','w')
    fp.write('# Standard deviation for variation over core (given as percentages of 1D values)\n')
    fp.write('# CoreRadiusFraction: {0:1.6f}\n'.format(coreRfraction))
    fp.write('#   Property\t')
    for case in DictOfCases.keys():
        fp.write('{0:s}\t'.format(case))
    fp.write('   max\t  mean\n')
    fp.write('#\n')
    for pr in exitVar:
        if pr not in ['dt_chem']:
            fp.write('{0:>12s}\t'.format(pr))
            x = np.array([])
            for case in DictOfCases.keys():
                x = np.append(x, std_as_percentage[case][pr])
                fp.write('{0:>6.4f}\t'.format(std_as_percentage[case][pr]))
            fp.write('{0:>6.4f}\t'.format(np.max(x)))
            fp.write('{0:>6.4f}\n'.format(np.mean(x)))
    fp.close()
    return 0

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr Examine Variation over core"
        print "    Version:", VERSION_STRING
        print "    To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
