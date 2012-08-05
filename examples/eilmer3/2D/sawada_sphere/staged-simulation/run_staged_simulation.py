#!/usr/bin/env python
# run_staged_simulation.py
#
# Top-level script to coordinate the running of the simulation in stages.
# The grid is refined with each stage.
#
# PJ, 13-Mar-2010 adapted from sphere-heat-transfer case
# RJG, 05-Aug-2012 adpated from lobb case

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_flow import StructuredGridFlow

# Allow script to access HOME variable globally
HOME = os.path.expandvars("$HOME")

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
    run_command("%s/e3bin/e3prep.py --job=%s --do-svg" % (HOME, stageName,))
    run_command("mpirun -np %d %s/e3bin/e3mpi.exe --job=%s --run" 
                % (np, HOME, stageName,))
    return

#---------------------------------------------------------------
def locate_shock_along_strip(x, y, p):
    """
    Shock location is identified as a pressure rise 
    along a strip of points.
    """
    n = len(x)
    p_max = max(p)
    p_trigger = p[0] + 0.3 * (p_max - p[0])
    x_old = x[0]; y_old = y[0]; p_old = p[0]
    for i in range(1,n):
        x_new = x[i]; y_new = y[i]; p_new = p[i]
        if p_new > p_trigger: break
        x_old = x_new; y_old = y_new; p_old = p_new
    frac = (p_trigger - p_old) / (p_new - p_old)
    x_loc = x_old * (1.0 - frac) + x_new * frac
    y_loc = y_old * (1.0 - frac) + y_new * frac
    return x_loc, y_loc

def locate_shock_front(stageName, nbi, nbj):
    """
    Reads all flow blocks and returns the coordinates 
    of the shock front, searching along the stagnation line.
    """
    blockData = []
    for ib in range(nbi):
        blockData.append([])
        for jb in range(nbj):
            blkindx = ib*nbj + jb
            fileName = 'flow/t9999/%s.flow.b%04d.t9999.gz' \
                % (stageName, blkindx)
            fp = gzip.open(fileName, "r")
            blk = StructuredGridFlow()
            blk.read(fp)
            blockData[ib].append(blk)
            fp.close()
    jb = 0
    nj = blockData[0][jb].nj
    j = 0
    x = []; y = []; p = [];
    for ib in range(nbi):
        ni = blockData[ib][jb].ni
        k = 0 # 2D only
        for i in range(ni):
            x.append(blockData[ib][jb].data['pos.x'][i,j,k])
            y.append(blockData[ib][jb].data['pos.y'][i,j,k])
            p.append(blockData[ib][jb].data['p'][i,j,k])
    return locate_shock_along_strip(x, y, p)

#---------------------------------------------------------------
# main...

jobName = 'ss3'
Rc = 31.8e-3 # m, radius of sphere
Db = 2.0*Rc

# Free-stream flow definition,
p_inf = 20.0e3  # Pa 
T_inf = 296.0   # K
u_inf = 4.68e3  # m/s

nbi = 2; nbj = 2
np = nbi * nbj  # number of processes for MPI simulation
# starting values for stage 0
paramDict = {'jobName': jobName, 'stage':0, 'Rc':Rc, 'flow_lengths':30,
             'p_inf':p_inf, 'T_inf':T_inf, 'u_inf':u_inf,
             'nn':60, 'nbi':nbi, 'nbj':nbj, 'np':np}

# We will work through the following resolutions, one for each stage.
nn_list = [30, 60, 90, 135]
fl_list = [30, 3, 2, 2]
# ...and collect the d/D values for a summary at the end.
d_list = []

for stage in range(len(nn_list)):
    paramDict['stage'] = stage
    paramDict['nn'] = nn_list[stage]
    paramDict['flow_lengths'] = fl_list[stage]
    run_stage(paramDict, jobName, stage)
    x_shock, y_shock = locate_shock_front(jobName+str(stage), nbi, nbj)
    d_list.append(-(x_shock+Rc))
    print "stage=", stage, "x_shock=", x_shock, "d=", d_list[-1]

print "Done at top-level."
print "Experimental result for shock detachment distance:"
print "    p=200.e3 Pa, T=296.0 K, u=4.68e3 m/s"
print "    d = 2.59 mm"
print " Sawada & Dendou CFD result: "
print "    d = 2.56 mm"
print "Simulation results:"
for stage in range(len(nn_list)):
    print "    nn=", nn_list[stage], "d=", d_list[stage]

