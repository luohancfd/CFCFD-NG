#! /usr/bin/env python 
# optimise_nozzle.py
#
# This script performs an optimisation of a given supersonic nozzle contour
# using the design methodology suggested by Chris Craddock (UQ Mechanical 
# Engineering Departmental Report 2/00). It couples the Nelder & Mead 
# optimiser (see the minimize function in cfcfd3/lib/cfpylib/nm/nelmin.py) 
# to a flow solver built specifically to run supersonic nozzle simulations, 
# NENZFr (see cfcfd3/app/nenzfr).
#
# In this script, the nozzle contour is optimised to produce a core flow
# at the exit of the nozzle that has the least deviations from the design 
# Mach number and has the least flow angularity. 
# 
# Wilson Chan, 31 Aug 2012.
# ----------------------------------------------------------------------------------


import sys, os
from math import tan, radians
import shlex, subprocess
from subprocess import PIPE
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from cfpylib.nm.nelmin import *


def run_command(cmdText):
    """
    Run the command as a subprocess.
    """
    print "About to run cmd:", cmdText
    if (type(cmdText) is list):
        args = cmdText
    else:
        args = shlex.split(cmdText)
    p = subprocess.Popen(args)
    # wait until the subprocess is finished
    stdoutData, stderrData = p.communicate()
    return


def objective_function(y):
    """
    Objective function for the optimisation of the internal contour of a T4 
    nozzle. The input to this function, "y", is a list of the relative radial 
    coordinates of the Bezier control points for the nozzle. Each radial 
    coordinate in this list is relative to the radial coordinate of the
    upstream point.
    """

    nFixedPts = 2  # Number of fixed Bezier control points

    # Target parameters
    M_target = 7.0        # Target Mach number
    dtheta_target = 0.02  # Target variation in outflow angle (in degrees)
    dM_target = 0.01      # Target variation in Mach number

    # Step 1: Read in the x- and y-coordinates of the Bezier control points of the 
    # initial nozzle contour.
    x_init = []; y_init = []
    fi = open('Bezier-control-pts-t4-m7.initial.data', 'r')
    fi.readline()  # Skip the first header line
    while True:
        buf = fi.readline().strip()
        if len(buf) == 0: break
        tokens = [float(word) for word in buf.split()]
        x_init.append(tokens[0])
        y_init.append(tokens[1])
    fi.close()
 
    # Compile a data file containing Bezier control points that will get
    # read in for a NENZFr simulation in the next step.
    fo = open('Bezier-control-pts-t4-m7.data', 'w')
    fo.write('#    x, m    y, m \n')
    # Write the coordinates for the control points that are fixed.
    for indx in range(nFixedPts):
        fo.write('%.7f %.7f \n' % (x_init[indx], y_init[indx]))
    # y_absolute is the y-coordinate of the first reference point.
    y_absolute = y_init[nFixedPts-1]
    # For the control points that have been shifted by the optimiser, we
    # use y_absolute to work out the absolute y-coordinate of each point
    # from the given relative coordinate in the "y" list.
    for indx in range(nFixedPts, len(y)+nFixedPts):
        y_absolute += y[indx-nFixedPts]
        fo.write('%.7f %.7f \n' % (x_init[indx], y_absolute))
    fo.close()

    # Run NENZFr in block-marching mode for the given nozzle contour.
    run_command('./run-nenzfr-in-mpi-block-marching-mode.sh')

    # Read in an extracted slice of data at the exit of the nozzle
    fi = open('nozzle-exit.data', 'r')
    # Keep a list of variables in order of appearance (copied from nenzfr.py).
    varLine = fi.readline().strip(); items = varLine.split()
    if items[0] == '#': del items[0]
    if items[0] == 'Variables:': del items[0]
    variable_list = [item.split(':')[1] for item in items]
    # Store the data in lists against these names (copied from nenzfr.py).
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

    # Compute secondary functions that contribute to the final objective function
    # (as per Equation 3 in Craddock's report).
    f_theta = 0.0; f_M = 0.0  # Initialise both functions to zero first.
    N = 0  # Counter for the number of cells in the core flow.
    for j in range(len(data['pos.y'])):
        # Compute the dMdy parameter. Note that dMdy is the radial gradient of Mach 
        # number, and is the parameter used by Chris Craddock to judge whether each
        # cell is in the core flow or not.
        if j == 0:
            # If we are in the first cell, set dMdy to a number larger than -20.
            # This is done to ensure that the first cell is always in the core flow.
            dMdy = 0.0  
        else:
            # Compute dMdy. 
            dMdy = (data['M_local'][j] - data['M_local'][j-1]) /\
                   (data['pos.y'][j] - data['pos.y'][j-1])
        if dMdy >= -20.0:
            # This cell is in the core flow, so we use the cell-local flow parameters 
            # (e.g. vel.x, vel.y and Mach number) to compute and to accumulate values 
            # for the secondary functions.
            f_theta += (data['vel.y'][j] / data['vel.x'][j])**2
            f_M += (data['M_local'][j] - M_target)**2
            N += 1

    # Weighting parameters (as per Equation 4 in Craddock's report).
    phi_theta = 1.0 / tan(radians(dtheta_target))
    phi_M = 1.0 / dM_target

    # Weight the secondary functions by weighting parameters (as per 
    # Equation 3 in Craddock's report).
    f_theta = phi_theta**2 / N * f_theta
    f_M = phi_M**2 / N * f_M

    # Objective function (as per Equation 2 in Craddock's report).
    obj_funct = (f_theta + f_M)**2

    return obj_funct


def get_design_variables(filename, nFixedPts):
    """
    Generate the design variables needed for the optimiser from a initial nozzle
    contour. The inputs to this function are
       1. the name of the file containing the Bezier control points of the 
          initial nozzle contour and
       2. the number of control points in this file that will not be varied by
          the optimiser.
    This function returns a list of y-coordinates of the Bezier control points
    that will be varied by the optimiser. Note that each of these y-coordinates
    are not absolute y-coordinates, but are instead y-coordinates that are 
    relative to the upstream Bezier control point (as specified in Section 3 of
    Craddock's report).
    """
    y_absolute = []; y_relative = []
    fi = open(filename, 'r')
    fi.readline()  # Skip the first header line.
    # Read and store all absolute y-coordinates.
    while True:
        buf = fi.readline().strip()
        if len(buf) == 0: break
        tokens = [float(word) for word in buf.split()]
        y_absolute.append(tokens[1])
    fi.close()
    # Compute and store the relative y-coordinates only for the control points
    # that will be varied by the optimiser.
    for indx in range(nFixedPts, len(y_absolute)):
        # The relative y-coordinate is computed by subtracting the y-coordinate
        # for the current control point by that of the previous control point.
        y_relative.append(y_absolute[indx] - y_absolute[indx-1])
    return y_relative


# ----------------------------------------------------------------------------------

if __name__ == '__main__':

    print "Begin the optimisation of the nozzle contour..."
    print "---------------------------------------------------"

    # Generate design variables from the initial given nozzle contour.
    x = get_design_variables('Bezier-control-pts-t4-m7.initial.data', 2)

    # Increments to apply to x when forming the initial simplex.
    dx = [0.003, 0.003, 0.005, 0.005, 0.005, 0.005, 0.005]

    # Check that the sizes for the x and dx lists are the same.
    if len(x) != len(dx):
        print "Error! x and dx lists must have the same number of variables."
        print "x has ", len(x), "variables, while dx has ", len(dx), "variables."
        sys.exit()

    # Input parameters for the optimiser.
    tol = 1.0     # Value at which optimisation convergence is achieved.
    maxfe = 1000  # Maximum number of function evaluations.
    n_check = 5   # Number of steps between convergence checks.

    # Start optimization.
    x, fx, conv_flag, nfe, nres = minimize(objective_function, x, dx, tol, maxfe, n_check)

    # Print the results and optimisation statistics that have been returned by the 
    # optimiser when it reaches convergence.
    print "x = ", x
    print "fx =", fx
    print "convergence-flag =", conv_flag
    print "number-of-fn-evaluations =", nfe
    print "number-of-restarts =", nres
    print "---------------------------------------------------"
    print "Done."
