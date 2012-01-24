#!/usr/bin/env python
# nenzfr.py
# NENZF reloaded (Non-Equilibrium Nozzle Flow reloaded)
#
# This script coordinates the running of the T4 nozzle calculation.
# The intention is to provide a fairly quick calculation of 
# the test flow conditions at the exit plane of the selected nozzle.
# Behind the scene, estcj.py is used to get an estimate of 
# the flow condition at the nozzle throat and then Eilmer3 is used
# to compute the expanding flow in the divergent part of the nozzle.
# Finally, a profile is examined at the downstream-end of the nozzle
# to extract nominal flow condition data.
# 
# Peter Jacobs, Luke Doherty, Wilson Chan and Rainer Kirchhartz
# School of Mechancial and Mining Engineering
# The University of Queensland

VERSION_STRING = "05-Oct-2011"

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
import optparse
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#---------------------------------------------------------------

def prepare_input_script(substituteDict, jobName):
    """
    Prepare the actual input file for Eilmer3 from a template.
    """
    templateFileName = E3BIN+"/nenzfr_data_files/nozzle.input.template"
    scriptFileName = jobName + ".py"
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
    # wait until the subprocess is finished
    stdoutData, stderrData = p.communicate() 
    return

def quote(str):
    """
    Put quotes around a string.
    """
    return '"' + str + '"'

def print_stats(sliceFileName,jobName):
    """
    Display statistics of flow properties at the nozzle exit.
    """
    print "Nozzle-exit statistics:"
    fp = open(sliceFileName, 'r')
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
    #
    # Identify edge of core flow.
    ys = data['pos.y']
    y_edge = ys[-1] * 2.0/3.0  # somewhat arbitrary...
    #
    # Compute and print area-weighted-average core flow values.
    exclude_list = ['pos.x', 'pos.y', 'pos.z', 'volume', 'vel.z', 'S']
    #
    fout = open(jobName+'-exit.stats','w')
    fout.write('%10s  %12s   %10s  %10s %10s\n' % \
                   ("variable","mean-value","minus","plus","std-dev"))
    fout.write(60*'-')
    fout.write('\n')
    #
    print "%10s  %12s    %10s %10s %10s" % \
        ("variable", "mean-value", "minus", "plus","std-dev")
    print 60*'-'
    for var in variable_list:
        if var in exclude_list: continue
        A = 0.0; F = 0.0;
        for j in range(len(ys)):
            if ys[j] > y_edge: break
            if j == 0:
                y0 = 0.0
            else:
                y0 = 0.5*(ys[j-1]+ys[j])
            y1 = 0.5*(ys[j]+ys[j+1])
            dA = y1**2 - y0**2
            F += data[var][j] * dA
            A += dA
        mean = F/A
        # Identify low and high values.
        diff_minus = 0.0
        diff_plus = 0.0
        count = 0.0
        stddev = 0.0
        for j in range(len(ys)):
            if ys[j] > y_edge: break
            diff = data[var][j] - mean
            diff_minus = min(diff, diff_minus)
            diff_plus = max(diff, diff_plus)
            count += 1
            stddev += diff**2
        # Calculate the sample standard deviation
        stddev = (stddev/(count-1))**0.5
        print "%10s  %12.4g    %10.3g %10.3g %10.3g" % \
              (var, mean, diff_minus, diff_plus, stddev)
        fout.write('%10s  %12.4g    %10.3g %10.3g %10.3g\n' % \
              (var, mean, diff_minus, diff_plus, stddev))
    #    
    print 60*'-'
    #
    fout.write(60*'-')
    fout.close()
    return


def main():
    """
    Examine the command-line options to decide the what to do
    and then coordinate the calculations done by Eilmer3.
    """
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('--gas', dest='gasName', default='air',
                  choices=['air', 'air5species', 'n2', 'co2', 'h2ne'],
                  help=("name of gas model: "
                        "air; " "air5species; " "n2; " "co2; " "h2ne"))
    op.add_option('--p1', dest='p1', type='float', default=None,
                  help=("shock tube fill pressure or static pressure, in Pa"))
    op.add_option('--T1', dest='T1', type='float', default=None,
                  help=("shock tube fill temperature, in degrees K"))
    op.add_option('--Vs', dest='Vs', type='float', default=None,
                  help=("incident shock speed, in m/s"))
    op.add_option('--pe', dest='pe', type='float', default=None,
                  help=("equilibrium pressure (after shock reflection), in Pa"))
    op.add_option('--chem', dest='chemModel', default='eq',
                  choices=['eq', 'neq', 'frz', 'frz2'], 
                  help=("chemistry model: " "eq=equilibrium; " 
                        "neq=non-equilibrium; " "frz=frozen " 
                        "[default: %default]"))
    op.add_option('--area', dest='areaRatio', default=27.0,
                  help=("nozzle area ratio. only used for estcj calc. "
                        "use when --cfile(--gfile) are "
                        "specified. [default: %default]"))
    op.add_option('--job', dest='jobName', default='nozzle',
                  help="base name for Eilmer3 files [default: %default]")
    op.add_option('--cfile', dest='contourFileName', 
                  default='contour-t4-m4.data',
                  help="file containing nozzle contour [default: %default]")
    op.add_option('--gfile', dest='gridFileName', default='None',
                  help="file containing nozzle grid. "
                  "overrides --cfile if both are given "
                  "[default: %default]")
    op.add_option('--exitfile', dest='exitSliceFileName', 
                  default='nozzle-exit.data',
                  help="file for holding the nozzle-exit data [default: %default]")
    op.add_option('--just-stats', dest='justStats', action='store_true', 
                  default=False,
                  help="skip the detailed calculations and "
                  "just retrieve exit-flow statistics")
    # The following defaults suit Like's Mach 10 calculations.
    op.add_option('--nni', dest='nni', type='int', default=1800,
                  help=("number of axial cells"))
    op.add_option('--nnj', dest='nnj', type='int', default=300,
                  help=("number of radial cells"))
    op.add_option('--nbi', dest='nbi', type='int', default=180,
                  help=("number of axial blocks for the divergence section (nozzle_blk)"))
    op.add_option('--bx', dest='bx', type='float', default=1.05,
                  help=("clustering in the axial direction"))
    op.add_option('--by', dest='by', type='float', default=1.005,
                  help=("clustering in the radial direction"))
    op.add_option('--max-time', dest='max_time', type='float', default=6.0e-3,
                  help=("overall simulation time for nozzle flow"))
    op.add_option('--max-step', dest='max_step', type='int', default=80000,
                  help=("maximum simulation steps allowed"))
    opt, args = op.parse_args()
    #
    # If we have already run a calculation, it may be that we just want
    # to extract the exit-flow statistics again.
    if opt.justStats:
        print_stats(opt.exitSliceFileName,opt.jobName)
        return 0
    #
    # Go ahead with a new calculation.
    # First, make sure that we have the needed parameters.
    bad_input = False
    if opt.p1 is None:
        print "Need to supply a float value for p1."
        bad_input = True
    if opt.T1 is None:
        print "Need to supply a float value for T1."
        bad_input = True
    if opt.Vs is None:
        print "Need to supply a float value for Vs."
        bad_input = True
    if opt.pe is None:
        print "Need to supply a float value for pe."
        bad_input = True
    if bad_input:
        return -2
    #
    # Get the nozzle contour file into the current work directory.
    run_command('cp '+E3BIN+'/nenzfr_data_files/'+opt.contourFileName+' .')
    # Set up the equilibrium gas-model file as a look-up table.
    eqGasModelFile = 'cea-lut-'+opt.gasName+'.lua.gz'
    if not os.path.exists(eqGasModelFile):
        run_command('build-cea-lut --case='+opt.gasName)
    # Runs estcj to get the equilibrium shock-tube conditions up to the nozzle-supply region.
    command_text = E3BIN+('/estcj.py --gas=%s --T1=%g --p1=%g --Vs=%g --pe=%g --task=st --ofn=%s' % 
                          (opt.gasName, opt.T1, opt.p1, opt.Vs, opt.pe, opt.jobName))
    run_command(command_text)
    if opt.chemModel in ['eq',]:
        # The gas model description for Eilmer3 is entirely in the look-up table file.
        gmodelFile = eqGasModelFile
    else:
        # We'll assume that the gas-model file of default name is set up.
        # TODO: Luke, this needs to be modified, I suspect.
        gmodelFile = 'gas-model.lua'
    # Set up the input script for Eilmer3.
    paramDict = {'jobName':quote(opt.jobName), 'gasName':quote(opt.gasName),
                 'T1':opt.T1, 'p1':opt.p1, 'Vs':opt.Vs, 'pe':opt.pe,
                 'contourFileName':quote(opt.contourFileName),
                 'gridFileName':quote(opt.gridFileName), 
                 'chemModel':quote(opt.chemModel),
                 'areaRatio':opt.areaRatio,
                 'nni':opt.nni, 'nnj':opt.nnj, 'nbi':opt.nbi, 'bx':opt.bx, 'by':opt.by,
                 'max_time':opt.max_time, 'max_step':opt.max_step}
    prepare_input_script(paramDict, opt.jobName)
    # Run Eilmer3
    run_command(E3BIN+('/e3prep.py --job=%s --do-svg' % (opt.jobName,)))
    run_command(E3BIN+('/e3shared.exe --job=%s --run' % (opt.jobName,)))
    # Exit plane slice
    run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
                +('--output-file=%s --slice-list="-1,-2,:,0" ' % 
                  (opt.exitSliceFileName,))
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Centerline slice
    run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
                +('--output-file=%s-centreline.data --slice-list=":,:,1,0" ' % 
                  (opt.jobName,))
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Prep files for plotting with Paraview            
    run_command(E3BIN+('/e3post.py --job=%s --vtk-xml --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Generate averaged exit flow properties
    print_stats(opt.exitSliceFileName,opt.jobName)                
    # Compute viscous data at the nozzle wall
    run_command(E3BIN+'/nenzfr_compute_viscous_data.py --job=%s' % (opt.jobName,))

    #
    # The following are additional commands specific to Luke D. and the Mach 10 nozzle.
    # TODO: Luke, we need to think of a good way to encode all of the options.
    #       I've added quite a few command-line options and more are needed.
    #       Maybe a dictionary of customised parameters for each case.
    #
    # Extract a slice from the last block along jk index directions at the i-index that
    # is closest to the point x=1.642, y=0.0, z=0.0
    run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
               +('--output-file=%s2 --slice-at-point="-1,jk,1.642,0,0" ' % 
                 (opt.exitSliceFileName,))
               +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')

    # 29-07-2011 Luke D. 
    # Copy the original exit stats file to a temporary file
    run_command(('cp %s-exit.stats %s-exit.stats_temp') % (opt.jobName, opt.jobName,))
    # (Re) Generate the exit stats using the slice-at-point data
    print_stats(opt.exitSliceFileName+'2',opt.jobName)
    run_command('cp %s-exit.stats %s-exit.stats2' % (opt.jobName, opt.jobName,))
    # Now rename the temporary exit stats file back to its original name
    run_command('mv %s-exit.stats_temp %s-exit.stats' % (opt.jobName, opt.jobName,))
    #
    # End specific commands
    # 
    return 0

#---------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr: Calculate Shock Tunnel Test Flow Conditions"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
