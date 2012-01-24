#! /usr/bin/env python
##
## \file design_x2_nozzle.py
## \ingroup geom
##
## \brief Use the sm_3d code in single-solution mode to do lots of
##        what-if simulations of a nozzle.
##
## \author PA Jacobs
##
## The nelmin optimizer supervises the calculations,
## with the goal of generating a good quality flow.
##
## \version 21-Mar-05 first attempt, inviscid design on Mach number only
##          22-Mar-05 inviscid design including both Mach number and angle
##

import os
from math import sqrt
from Numeric import array
from wall_definition import bezier_nodes_for_wall

def prepare_input_file(jobName, opt_parameters):
    """Prepares the .par file frrom a template and the supplied values.

    opt_parameters is the list of optimization parameters.
    See below to find their interpretation.
    """
    template = """The parameter file for ideal-gas design of X2 M10 contour.
0                            case_id - generic case
0.4   300  0.0001            CFL, max_t_steps, closing tolerance
2 3 0.9                      Xorder, i_supress, p_safety
0.0  0.0001  1.0  10000      Xi_0, dXi, X_max, max_x_steps 2000
0.005                        dXi_plot
1                            slice_ident
40  2                        nny, nnz
0 10                         smooth_grid, smooth_iter
5  3  3  3                   bc_N, E, S, W
296.0 296.0 296.0 296.0      Twall_N, E, S, W
<rho> <ux> 0.0 0.0 <e> free stream rho, ux, uy, uz, e
0 1 1                        use_B_spline, bezier_box, 2d format
"""
    # Inflow properties
    M = 7.0      # Mach number
    T = 3000.0   # degrees K
    p = 100.0e3  # Pascals
    g = 1.4      # ratio of specific heats
    Rgas = 287.0 # gas constant J/(kg.K)
    Cv = Rgas / (g - 1.0)
    rho = p / (Rgas * T)
    e = Cv * T
    a = sqrt(g * Rgas * T)
    ux = M * a
    template = template.replace("<rho>", ("%12.4e" % rho))
    template = template.replace("<ux>", ("%12.4e" % ux))
    template = template.replace("<e>", ("%12.4e" % e))
    #
    # Wall profile
    x_list, y_list = bezier_nodes_for_wall(opt_parameters)
    n_nodes = len(x_list)
    template += ("%d\n" % n_nodes)
    for i in range(n_nodes):
        template += ("%12.4e %12.4e   A %d\n%12.4e %12.4e   B %d\n" % \
                     (x_list[i], 0.0, i, x_list[i], y_list[i], i))
    f = open(jobName+".par", "wt")
    f.write(template)
    f.close()
    return

def run_sm_prep(jobName, verbose=0):
    """Runs the inflow profile preparation program"""
    if verbose:
        print "Prepare inflow data."
    fi, fo = os.popen2("sm_prep.x")
    fi.write(jobName + "\n")
    resultFlag = fi.close()
    if resultFlag:
        print "Some problem with sm_prep:", resultFlag
    prep_stdout = fo.read()
    if verbose:
        print prep_stdout
    fo.close()
    return

def run_sm_3d_simulation(jobName, verbose=0):
    """Runs the sm_3d simulation program"""
    if verbose:
        print "Run simulation..."
    fi, fo = os.popen2("sm_3d.x")
    fi.write(jobName + "\n")
    resultFlag = fi.close()
    if resultFlag:
        print "Some problem with sm_3d:", resultFlag
    simulation_stdout = fo.read()
    if verbose:
        print simulation_stdout
    fo.close()
    return

def extract_profile_data(jobName, verbose=0):
    """Runs the postprocessing program.

    This should result in a .plt file being written into the
    current directory.
    The TECPLOT format of the first few lines is...

    TITLE = "Solution data from test " 
    VARIABLES =  "y"  ,"rho"  ,"p"  ,"M"  ,"T"  ,"Ang(deg)" 
    ZONE T="Variable iy, ix = 200, iz = 1 ", I=40, F=POINT 
    3.757928e-03 0.0329487 18621.4 8.79282 1969.21 -0.930024 
    1.123198e-02 0.0231717 11191.2 9.55019 1682.83 -1.3846 
    """
    if verbose: print "Extracting profile...."
    fi, fo = os.popen2("sm_prof.x")
    fi.write(jobName +
"""
201 40 2
1
1 200 0 1
0 1 0
1 0 0 0 1 1 1 0 1
""")
    resultFlag = fi.close()
    if resultFlag:
        print "Some problem with sm_prof:", resultFlag
    prof_stdout = fo.read()
    if verbose:
        print prof_stdout
    fo.close()
    #
    # Now pick up the profile data.
    f = open(jobName + ".plt", "r")
    content = f.readlines()
    values = []
    for line in content[3:]:
        numbers = [float(word) for word in line.split()]
        values.append(numbers)
    va = array(values)
    y = va[:,0]
    rho = va[:,1]
    p = va[:,2]
    M = va[:,3]
    T = va[:,4]
    theta = va[:,5]
    return y, rho, p, M, T, theta

def objective_function(parameter_list, verbose=0):
    """Runs the simulation with current guess for parameters and
    returns an error indicator"""
    #
    if verbose: print "Begin objective function..."
    print "p_list=", [("%13.5e" % p) for p in parameter_list]
    jobName = "test"
    prepare_input_file(jobName, parameter_list)
    run_sm_prep(jobName, verbose)
    run_sm_3d_simulation(jobName, verbose)
    if verbose: print "End Sim Job"
    y, rho, p, M, T, theta = extract_profile_data(jobName, verbose)
    # print "y=", y, "M=", M
    #
    # Construct the objective function estimate from
    # (a) the deviation of mach number from the design value
    # (b) the deviation of flow angle from zero
    # across the exit plane.
    # For inviscid design, include the whole profile.
    # For viscous design, include only that part in the core flow.
    N = len(y)
    #
    M_design = 10.0   # desired Mach number
    dM_design = 0.01  # desired deviation limit
    dev = M - M_design
    f_M = 1.0/N * sum(dev * dev) / (dM_design**2)
    #
    theta_design = 0
    dtheta_design = 0.01
    dev = theta - theta_design
    f_theta = 1.0/N * sum(dev * dev) / (dtheta_design**2)
    #
    f_obj = (f_M + f_theta)**2
    #
    print "f_obj= %13.5e" % f_obj
    sys.stdout.flush()
    return f_obj

import nelmin
import sys

if __name__ == '__main__':
    # We are running a stand-alone script, so get on with some work.
    if len(sys.argv) == 1 or (len(sys.argv) > 1 and sys.argv[1]) == '-help':
        print "Usage: design_x2_nozzle.py [-opt|-single|-help]"
        sys.exit()

    print "Begin"
    # param = [0.06, 0.07, 0.08, 0.09, 0.10]
    param = [0.040496, 0.065839, 0.051086, 0.089875, 0.103871]
    if sys.argv[1] == '-opt':
        dparam = [0.002,] * len(param)  # nominal perturbations
        popt, fpopt, conv_flag, nfe, nres = \
           nelmin.minimize(objective_function, param, dparam, 1.0e-6, 100)
        print "optimised parameters=", popt
        print "objective=", fpopt
        print "convergence-flag=", conv_flag
        print "number-of-fn-evaluations=", nfe
        print "number-of-restarts=", nres
    elif sys.argv[1] == '-single':
        objective_function(param, 1)
    else:
        print "Unknown command option:", sys.argv[1]
    print "Done"
