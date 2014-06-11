#! /usr/bin/env python 
# optimise_T4_Mach_4_nozzle.py
#
# This script optimises an initial contour for a supersonic nozzle 
# using the design methodology suggested by Chris Craddock (UQ 
# Mechanical Engineering Departmental Report 2/00). It uses the
# minimize() function in cfcfd3/lib/cfpylib/nm/nelmin.py. 
#
# Wilson Chan, 30 Nov 2013.
#
# Wilson Chan, 05 Jun 2014 
#   - added a secondary function, as recommended by Han Wei, to stop
#     the nozzle contour from curving towards the nozzle axis. The 
#     inward-curving phenomena is brought about by the optimiser trying 
#     to create an oblique shock that cancels out the expansion waves
#     to get the coreflow more uniform. See f_penalty in the code.
# ------------------------------------------------------------------

import sys, os
from math import tan, radians
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from cfpylib.nm.nelmin import *
from nenzfr_utils import run_command

def objective_function(y):
    """
    Objective function for the design of the internal contour of 
    a T4 nozzle. The input y is a list of coordinates of the 
    Bezier control points for the nozzle. The radial coordinates 
    are specified relative to the radial coordinate of the 
    upstream point.
    """

    nFixedPts = 2  # Number of fixed Bezier control points

    # Targets
    M_target = 4.0        # Target Mach number
    dtheta_target = 0.02  # Target variation in outflow angle (in degrees)
    dM_target = 0.01      # Target variation in Mach number

    # Flag for the inclusion of a secondary penalty function in our optimisation. 
    # Switch this on, if your optimised nozzle contour always ends with a negative
    # gradient (i.e. the nozzle curves towards its axis near the exit of the contour.) 
    # This secondary penalty function should limit this inward-turning phenomena.
    include_penalty_function = 0

    # Read in the x-coordinates for all points and the y-coordinates for
    # the first nFixedPts points - they are fixed and not changed in the 
    # optimisation run.
    x_fixed = []; y_fixed = []
    fi = open(basename+'.initial.data', 'r')
    fi.readline()  # Skip the first line
    while True:
        buf = fi.readline().strip()
        if len(buf) == 0: break
        tokens = [float(word) for word in buf.split()]
        x_fixed.append(tokens[0])
        y_fixed.append(tokens[1])
    fi.close()
 
    # Use x and y to create new data file containing Bezier control
    # points to generate the internal contour of the nozzle.
    fo = open(basename+'.data', 'w')
    fo.write('#    x, m    y, m \n')
    for indx in range(nFixedPts):
        fo.write('%.7f %.7f \n' % (x_fixed[indx], y_fixed[indx]))
    y_absolute = y_fixed[nFixedPts-1]
    for indx in range(nFixedPts, len(y)+nFixedPts):
        y_absolute += y[indx-nFixedPts]
        fo.write('%.7f %.7f \n' % (x_fixed[indx], y_absolute))
    fo.close()

    # Run nenzfr with the given nozzle contour.
    run_command('./run-nenzfr-mpi-block-marching-mode.sh')

    # Read in an extracted slice of data at the exit of the nozzle
    fi = open('nozzle-exit.data', 'r')
    # Keep a list of variables in order of appearance (from nenzfr.py).
    varLine = fi.readline().strip(); items = varLine.split()
    if items[0] == '#': del items[0]
    if items[0] == 'Variables:': del items[0]
    variable_list = [item.split(':')[1] for item in items]
    # Store the data in lists against these names (from nenzfr.py).
    data = {}
    for var in variable_list:
        data[var] = []
    for line in fi.readlines():
        items = line.strip().split()
        if items[0] == '#': continue
        assert len(items) == len(variable_list)
        for i in range(len(items)):
            data[variable_list[i]].append(float(items[i]))    
    fi.close()

    # Secondary functions that contribute to the objective function.
    f_theta = 0.0; f_M = 0.0 # Initialise both functions to zero first.
    N = 0 # Initialise the counter for the number of cells in the boundary layer.
    for j in range(len(data['pos.y'])):
        # Definition used by Chris Craddock to estimate the boundary layer edge
        if j == 0: 
            dMdy = 0.0  # Set to some number so that the first point 
                        # is not set as the boundary layer edge.
        else:
            dMdy = (data['M_local'][j] - data['M_local'][j-1]) /\
                   (data['pos.y'][j] - data['pos.y'][j-1])
        # If dMdy >= -20.0, then we are in the core flow.
        if dMdy >= -20.0:
            f_theta += (data['vel.y'][j] / data['vel.x'][j])**2
            f_M += (data['M_local'][j] - M_target)**2
            N += 1

    # Weighting parameters.
    phi_theta = 1.0 / tan(radians(dtheta_target))
    phi_M = 1.0 / dM_target

    # Weight the secondary functions by weighting parameters.
    f_theta = phi_theta**2 / N * f_theta
    f_M = phi_M**2 / N * f_M

    # Secondary penalty function
    if include_penalty_function == 1:
        # Whenever the Bezier control points start going towards the nozzle axis, 
        # impose a really large penalty. Note though that this means that the 
        # optimiser might never reach its target objective of 1.0.
        if min(y) < 0.0:
            f_penalty = 1e9
        else:
            f_penalty = 0.0

    # Objective function
    if include_penalty_function == 1:
        obj_funct = (f_theta + f_M + f_penalty)**2
    else:
        obj_funct = (f_theta + f_M)**2

    return obj_funct


def get_design_variables(filename, nFixedPts):
    """
    Generate the design variables from a given file. The other 
    input is the number of fixed points in the data file that 
    are not part of the list of design variables for the
    optimisation run.

    Returns a list of y-coordinates of the Bezier control points
    that define the internal contour of a nozzle. Note that these
    coordinates are relative to those of the upstream point.
    """
    y_absolute = []; y_relative = []
    fi = open(filename, 'r')
    fi.readline()  # Skip the first line
    while True:
        buf = fi.readline().strip()
        if len(buf) == 0: break
        tokens = [float(word) for word in buf.split()]
        y_absolute.append(tokens[1])
    fi.close()
    for indx in range(nFixedPts, len(y_absolute)):
        y_relative.append(y_absolute[indx] - y_absolute[indx-1])
    return y_relative


# ------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin the optimisation of the T4 Mach 4 nozzle..."
    print "---------------------------------------------------"
    # Basename for file containing Bezier control points
    global basename
    basename = 'Bezier-control-pts-t4-m4'
    # List of N design variables 
    x = get_design_variables(basename+'.initial.data', 2)
    # List of N increments to apply to x when 
    # forming the initial simplex.
    dx = [0.003, 0.003, 0.005, 0.005, 0.005, 0.005, 0.005]
    # Check that the sizes for the x and dx lists are the same.
    if len(x) != len(dx):
        print "Error! x and dx lists must have the same number of variables."
        print "x has ", len(x), "variables, while dx has ", len(dx), "variables."
        sys.exit()
    # The terminating limit for the standard-deviation
    # the simplex function values.
    tol = 1.0
    # Maximum number of function evaluations 
    maxfe = 1000
    # Number of steps between convergence checks
    n_check = 5
    # Start function minimisation
    x, fx, conv_flag, nfe, nres = minimize(objective_function, x,
                                           dx, tol, maxfe, n_check)
    print "x = ", x
    print "fx =", fx
    print "convergence-flag =", conv_flag
    print "number-of-fn-evaluations =", nfe
    print "number-of-restarts =", nres
    print "---------------------------------------------------"
    print "Done."
