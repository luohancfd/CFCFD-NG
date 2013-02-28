#!/usr/bin/env python
# run_adaptive_simulation.py
#
# Top-level script to coordinate the running of the 
# shock-fitting simulation in stages.
#
# PJ, 22-Feb-2010
# AP, 17-Jan-2013 - Modified for shock-fitting.

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_flow import StructuredGridFlow
from math import sqrt, sin, cos, pi


#---------------------------------------------------------------
def prepare_input_script(substituteDict, jobName, stage):
    """
    Prepare the actual input file from a template.
    """
    stageName = jobName + str(stage)
    templateFileName = jobName + ".input.template"
    scriptFileName = stageName + ".py"
    fp = open(templateFileName, 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.substitute(substituteDict)
    fp = open(scriptFileName, 'w')
    fp.write(text)
    fp.close()
    return

def run_command(cmdText):
    """
    Run the command as a subprocess.
    """
    print "About to run cmd:", cmdText
    args = shlex.split(cmdText)
    p = subprocess.Popen(args)
    stdoutData, stderrData = p.communicate() 
    # wait until the subprocess is finished
    return

def run_stage(paramDict, jobName, stage):
    """
    Set up and run one stage of the simulation as a normal job.
    """
    prepare_input_script(paramDict, jobName, stage)
    stageName = jobName+str(stage)
    homePath = os.path.expanduser("~")
    run_command("%s/e3bin/e3prep.py --job=%s" % (homePath, stageName,))
    run_command("mpirun -np %d %s/e3bin/e3mpi.exe --job=%s --run" 
                % (np, homePath, stageName,))
    return

#---------------------------------------------------------------
# main...

jobName = 'sphere'

# Free-stream flow definition,
# We have initially static gas, processed by a Mach 5 shock.
# Inflow conditions are thus post-shock conditions.

M_inf = 5.0
R_gas = 287.1
g_gas = 1.4

T0_inf_R = 816.5                               # Freestream stagnation temperature, Rankine
T0_inf = T0_inf_R * 5.0/9.0                     # Freestream stagnation temperature, Kelvin
T_inf = T0_inf / ( 1 + 0.5 * (g_gas-1) * M_inf**2 )     # Freestream temperature, Kelvin
T_body = 0.218 * T0_inf                          # Body temperature, Kelvin

rho_inf = 1.71e4 * 515.378818     # (slug/ft^3 to kg/m^3) # density, kg/m^3
p_inf = 0.28 * 6894.75729       # (psi to Pa)           # freestream pressure, Pa
u_inf = M_inf * (g_gas * R_gas * T_inf)**0.5                    # flow speed, m/s

p_init = 0.3 * p_inf                                    # initial pressure, Pa

R = 2.5 * 0.0254 # (in to m)    # Sphere radius, m

# Initial simulation using guessed inflow boundary position.
stage = 0
x_d = [-1.5*R, -1.5*R, -1.0*R, 0.0]
y_d = [0.0, 1.0*R, 2.0*R, 3.0*R]
factor = 1; ni_basic = 12; nj_basic = 24 # Coarse grid to start to deal with instability.
nbi = 1; nbj = 4
np = nbi * nbj  # number of processes for MPI simulation
np0 = np
paramDict = {'jobName': jobName, 'stage':stage, 
             'R':R, 'x_d':x_d, 'y_d':y_d,
             'T_body':T_body, 'p_init':p_init,
             'p_inf':p_inf, 'T_inf':T_inf, 'u_inf':u_inf,
             'ni':factor*ni_basic, 'nj':factor*nj_basic,
             'nbi':nbi, 'nbj':nbj, 'np':np, 'np0':np0,
             'viscous_flag':0, 'viscous_delay':0.0,
             'shock_fitting_bc':False,'clustering':False,
             'body_lengths':20,'use_same_grid':True,'shock_fitting_flag':1}

run_stage(paramDict, jobName, stage)

# Start moving shock.
paramDict['shock_fitting_bc'] = True

stage += 1
paramDict['stage'] = stage
factor = 1
paramDict['ni'] = factor*ni_basic
paramDict['nj'] = factor*nj_basic
paramDict['body_lengths'] = 20.0
run_stage(paramDict, jobName, stage)
print "Stage ", stage

# Turn on viscosity and clustering, and refine grid,
# starting from previous solution.
paramDict['viscous_flag'] = 1
paramDict['clustering'] = True

stage += 1
paramDict['stage'] = stage
factor = 1
paramDict['ni'] = factor*ni_basic
paramDict['nj'] = factor*nj_basic
paramDict['body_lengths'] = 10.0
run_stage(paramDict, jobName, stage)
print "Stage ", stage

# Refine grid.
stage += 1
paramDict['stage'] = stage
factor = 2
paramDict['ni'] = factor*ni_basic
paramDict['nj'] = factor*nj_basic
paramDict['body_lengths'] = 10.0
run_stage(paramDict, jobName, stage)
print "Stage ", stage


