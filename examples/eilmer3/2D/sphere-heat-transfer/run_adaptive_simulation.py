#!/usr/bin/env python
# run_adaptive_simulation.py
#
# Top-level script to coordinate the running of the 
# solution-adaptive simulation in stages.
# This approximates a form of solution adaptivity in that
# the grid is adjusted to the shock occasionally.
# The grid is also refined with the stages.
#
# PJ, 22-Feb-2010

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_flow import StructuredGridFlow

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
    run_command("/home/peterj/e3bin/e3prep.py --job=%s --do-svg" % (stageName,))
    run_command("mpirun -np %d /home/peterj/e3bin/e3mpi.exe --job=%s -q --run" 
                % (np, stageName,))
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
    of the shock front in lists of coordinates.
    """
    blockData = []
    for ib in range(nbi):
        blockData.append([])
        for jb in range(nbj):
            blkindx = ib*nbj + jb
            fileName = 'flow/t0005/%s.flow.b%04d.t0005.gz' \
                % (stageName, blkindx)
            fp = gzip.open(fileName, "r")
            blockData[ib].append(StructuredGridFlow())
            blockData[ib][-1].read(fp)
            fp.close()
    x_shock = []; y_shock = []
    for jb in range(nbj):
        nj = blockData[0][jb].nj
        for j in range(nj):
            x = []; y = []; p = [];
            for ib in range(nbi):
                ni = blockData[ib][jb].ni
                k = 0 # 2D only
                for i in range(ni):
                    x.append(blockData[ib][jb].data['pos.x'][i,j,k])
                    y.append(blockData[ib][jb].data['pos.y'][i,j,k])
                    p.append(blockData[ib][jb].data['p'][i,j,k])
            xshock, yshock = locate_shock_along_strip(x, y, p)
            x_shock.append(xshock)
            y_shock.append(yshock)
    return x_shock, y_shock

#---------------------------------------------------------------
def define_bezier_points(alpha, x_s, y_s):
    """
    It is assumed that the centre of the circular body is at (0,0)
    and that we have a third-order Bezier curve that goes through
    the start and finish of the shock.
    """
    import math
    x0 = x_s[0]; y0 = 0.0         # the first point coincides with the shock
    x3 = 0.0; y3 = y_s[-1]        # locate final point also on shock
    x1 = x0; y1 = 0.5 * y3
    L = 0.4 * y3
    x2 = x3 - L * math.cos(alpha)
    y2 = y3 - L * math.sin(alpha)
    return [x0, x1, x2, x3], [y0, y1, y2, y3]

def fit_bezier(x_s, y_s):
    """
    Fits a Bezier curve to the shock coordinates 
    and returns lists of coordinates.
    """
    from cfpylib.nm.line_search import minimize
    import math
    #
    def objective(alpha, x_s=x_s, y_s=y_s):
        """
        Objective function for the optimizer.
        """
        from libprep3 import Bezier, Vector
        bx, by = define_bezier_points(alpha, x_s, y_s)
        bpath = Bezier([Vector(bx[0],by[0]),Vector(bx[1],by[1]),
                        Vector(bx[2],by[2]),Vector(bx[3],by[3])])
        nbez = 1000
        pbez = []
        for i in range(nbez):
            t = 1.0/nbez * i
            pbez.append(bpath.eval(t))
        n = len(x_s)
        sum_sq_err = 0.0
        for j in range(n):
            min_dist = (x_s[j]-pbez[0].x)**2 + (y_s[j]-pbez[0].y)**2
            for i in range(1,nbez):
                dist = (x_s[j]-pbez[i].x)**2 + (y_s[j]-pbez[i].y)**2
                if dist < min_dist: min_dist = dist
            sum_sq_err += min_dist
        # print "alpha=", alpha, "sum_sq_err=", sum_sq_err
        return sum_sq_err
    #
    alphaL, alphaR = minimize(objective, 0.0, math.pi/4)
    best_alpha = 0.5*(alphaL+alphaR)
    return define_bezier_points(best_alpha, x_s, y_s)

#---------------------------------------------------------------
# main...

jobName = 'sphere'
R = 6.6e-3                   # nose radius of sphere
T_body = 296.0               # surface T

# Free-stream flow definition,
# We have initially static gas, processed by a Mach 8 shock.
# Inflow conditions are thus post-shock conditions.
p_init = 6.7                 # Pa
p_inf = 535.6                # Pa
T_inf = 2573.5               # degrees K
u_inf = 2436.5               # flow speed, m/s

# Initial simulation using guessed inflow boundary position.
stage = 0
x_d = [-1.5*R, -1.5*R, -1.0*R, 0.0]
y_d = [0.0, 1.0*R, 2.0*R, 3.0*R]
factor = 2; ni_basic = 10; nj_basic = 10
nbi = 2; nbj = 2
np = nbi * nbj  # number of processes for MPI simulation
paramDict = {'jobName': jobName, 'stage':stage, 
             'R':R, 'x_d':x_d, 'y_d':y_d,
             'T_body':T_body, 'p_init':p_init,
             'p_inf':p_inf, 'T_inf':T_inf, 'u_inf':u_inf,
             'ni':factor*ni_basic, 'nj':factor*nj_basic,
             'nbi':nbi, 'nbj':nbj, 'np':np,
             'viscous_flag':1, 'viscous_delay':10.0e-6, 
             'body_lengths':10}  # 20 body_lengths normally
run_stage(paramDict, jobName, stage)

# Restart from stage 0 flow data,
# bringing grid in closer to the shock.
stage = 1
paramDict['stage'] = stage
x_shock, y_shock = locate_shock_front(jobName+str(stage-1), nbi, nbj)
x_d, y_d = fit_bezier(x_shock, y_shock)
# Scale out so that we are sure to capture the shock.
paramDict['x_d'] = [1.05*x for x in x_d]
paramDict['y_d'] = [1.1*y for y in y_d]
paramDict['viscous_delay'] = 0.0
run_stage(paramDict, jobName, stage)

# Restart from stage 1 flow data, refining grid
stage = 2
paramDict['stage'] = stage
factor = 3
paramDict['ni'] = factor*ni_basic
paramDict['nj'] = factor*nj_basic
paramDict['body_lengths'] = 5.0
run_stage(paramDict, jobName, stage)

# Restart from stage 2 flow data, 
# adjusting the inflow boundary and refining grid
stage = 3
paramDict['stage'] = stage
x_shock, y_shock = locate_shock_front(jobName+str(stage-1), nbi, nbj)
x_d, y_d = fit_bezier(x_shock, y_shock)
# Scale out so that we are sure to capture the shock.
paramDict['x_d'] = [1.05*x for x in x_d]
paramDict['y_d'] = [1.1*y for y in y_d]
factor = 4
paramDict['ni'] = factor*ni_basic
paramDict['nj'] = factor*nj_basic
paramDict['body_lengths'] = 5.0
run_stage(paramDict, jobName, stage)

# Restart from stage 3 flow data, refining the grid only.
stage = 4
paramDict['stage'] = stage
factor = 8
paramDict['ni'] = factor*ni_basic
paramDict['nj'] = factor*nj_basic
paramDict['body_lengths'] = 5.0
run_stage(paramDict, jobName, stage)

print "Done at top-level."

