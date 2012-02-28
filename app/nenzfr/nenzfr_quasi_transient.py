#!/usr/bin/env python
# nenzfr_quasi_transient.py
#
# This script coordinates the running of a series of nenzfr 
# calculations based on an input supply pressure trace, thus
# producing a 'quasi-transient' estimate of the variation of 
# the nozzle exit properties. 
# 
# Luke Doherty
# School of Mechancial and Mining Engineering
# The University of Queensland

VERSION_STRING = "24-Feb-2012"

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
import optparse
#import numpy
from numpy import array, mean, logical_and, zeros
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#---------------------------------------------------------------

#def prepare_input_script(substituteDict, jobName):
#    """
#    Prepare the actual input file for Eilmer3 from a template.
#    """
#    templateFileName = E3BIN+"/nenzfr_data_files/nozzle.input.template"
#    scriptFileName = jobName + ".py"
#    fp = open(templateFileName, 'r')
#    text = fp.read()
#    fp.close()
#    template = string.Template(text)
#    text = template.substitute(substituteDict)
#    fp = open(scriptFileName, 'w')
#    fp.write(text)
#    fp.close()
#    return

def prepare_run_script(substituteDict, jobName):
    """
    Prepare the actual run file for Nenzfr from a template.
    """
    #templateFileName = E3BIN+"nenzfr_data_files/run_template.sh"
    templateFileName = "./run_template.sh"
    scriptFileName = "run_" + jobName + ".sh"
    fp = open(templateFileName, 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.substitute(substituteDict)
    fp = open(scriptFileName, 'w')
    fp.write(text)
    fp.close()
    return scriptFileName
    

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

def calculate_pe_values(pefile, tstart, nsteps, dt):
    """
    Calculate mean pe values for 'nsteps' intervals 
    of 'dt' starting at 'tstart' using the data in 
    "pefile". 
    """
    
    # Load data
    fp = open(pefile, 'r')
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
    
    # Put the data into numpy arrays
    time = array(data['time'])
    pressure = array(data['pe'])
    
    # Store the mean equilibrium pressures for each interval in an array
    pe = zeros(nsteps)
    for k in range(nsteps):
        pe[k] = mean(pressure[logical_and(time>=tstart+k*dt, time<=tstart+(k+1)*dt)])
        print "k= %d; tstart= %f; pe[k]= %f" % (k, tstart+k*dt, pe[k])
    
    return pe


def main():
    """
    Examine the command-line options to decide the what to do
    and then coordinate the calculations done by Eilmer3.
    """
    op = optparse.OptionParser(version=VERSION_STRING)

    op.add_option('--runCMD', dest='runCMD', default='qsub ',
                  help=("command used to execute the run-script file "
                        "[default: %default]"))
    op.add_option('--tstart', dest='tstart', type='float', default='1.5e-3',
                  help=("blah [default: %default]"))
    op.add_option('--nsteps', dest='nsteps', type='float', default=5,
                  help=("blah [default: %default]"))
    op.add_option('--dt', dest='dt', type='float', default=0.5e-3,
                  help=("blah [default: %default]"))

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
    
    #op.add_option('--pe', dest='pe', type='float', default=None,
    #              help=("equilibrium pressure (after shock reflection), in Pa"))
    op.add_option('--pefile', dest='peFileName', default=None,
                  help="file specifying the transient equilibrium pressure (in Pa) "
                       "[default: %default]")
    
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
                  help=("number of axial cells [default: %default]"))
    op.add_option('--nnj', dest='nnj', type='int', default=300,
                  help=("number of radial cells [default: %default]"))
    op.add_option('--nbi', dest='nbi', type='int', default=180,
                  help=("number of axial blocks for the divergence section (nozzle_blk) "
                  "[default: %default]"))
    op.add_option('--bx', dest='bx', type='float', default=1.05,
                  help=("clustering in the axial direction [default: %default]"))
    op.add_option('--by', dest='by', type='float', default=1.005,
                  help=("clustering in the radial direction [default: %default]"))
    op.add_option('--max-time', dest='max_time', type='float', default=6.0e-3,
                  help=("overall simulation time for nozzle flow [default: %default]"))
    op.add_option('--max-step', dest='max_step', type='int', default=80000,
                  help=("maximum simulation steps allowed [default: %default]"))
    opt, args = op.parse_args()
    #
    # If we have already run a calculation, it may be that we just want
    # to extract the exit-flow statistics again.
    #if opt.justStats:
    #    print_stats(opt.exitSliceFileName,opt.jobName)i
    #    return 
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
    if opt.peFileName is None:
        print "Need to supply a string value for pefile."
        bad_input = True
    if bad_input:
        return -2
    
    
    # Calculate the averaged equilibrium pressure for each
    # time interval...
    pe = calculate_pe_values(opt.peFileName, opt.tstart, opt.nsteps, opt.dt)
    
    #print "pe= ", pe
    #print "length(pe)= ", len(pe)
    #print "pe(2)= ", pe[2]
    
    # Set up the run script for Nenzfr.
    paramDict = {'jobName':quote(opt.jobName), 'gasName':quote(opt.gasName),
                 'T1':opt.T1, 'p1':opt.p1, 'Vs':opt.Vs,
                 'contourFileName':quote(opt.contourFileName),
                 'gridFileName':quote(opt.gridFileName),
                 'chemModel':quote(opt.chemModel),
                 'areaRatio':opt.areaRatio,
                 'nni':opt.nni, 'nnj':opt.nnj, 'nbi':opt.nbi, 'bx':opt.bx, 'by':opt.by,
                 'max_time':opt.max_time, 'max_step':opt.max_step,
                 'exitSliceFileName':opt.exitSliceFileName}
    
    #print "paramDict= ", paramDict
    #paramDictTemp = paramDict
    #paramDictTemp['pe'] = pe[2]
    #print "paramDictTemp= ", paramDictTemp

    for k in range(len(pe)):
        # Create sub-directory for the current case.
        run_command('mkdir ./case'+str(k))
        
        # Set up the run script for Nenzfr
        paramDict['pe'] = pe[k]
        scriptFileName = prepare_run_script(paramDict, opt.jobName+'_case'+str(k))
        
        # Move the run script to its sub-directory
        command_text = 'mv '+scriptFileName+' ./case'+str(k)+'/'+scriptFileName
        run_command(command_text)
        
        # Change into the sub-directory, ensure the run script is exectuable and
        # then run it
        os.chdir('case'+str(k))
        run_command('chmod u+x '+scriptFileName)
        os.system(opt.runCMD+scriptFileName)
        os.chdir('../')
    
    
    
    
    #
    # Get the nozzle contour file into the current work directory.
    #run_command('cp '+E3BIN+'/nenzfr_data_files/'+opt.contourFileName+' .')
    # Set up the equilibrium gas-model file as a look-up table.
    #eqGasModelFile = 'cea-lut-'+opt.gasName+'.lua.gz'
    #if not os.path.exists(eqGasModelFile):
    #    run_command('build-cea-lut --case='+opt.gasName)
    #
    #
    # Runs estcj to get the equilibrium shock-tube conditions up to the nozzle-supply region.
    #command_text = E3BIN+('/estcj.py --gas=%s --T1=%g --p1=%g --Vs=%g --pe=%g --task=st --ofn=%s' % 
    #                      (opt.gasName, opt.T1, opt.p1, opt.Vs, opt.pe, opt.jobName))
    #run_command(command_text)
    #
    #if opt.chemModel in ['eq',]:
    #    # The gas model description for Eilmer3 is entirely in the look-up table file.
    #    gmodelFile = eqGasModelFile
    #else:
    #    # We'll assume that the gas-model file of default name is set up.
    #    # TODO: Luke, this needs to be modified, I suspect.
    #    gmodelFile = 'gas-model.lua'
    


    # Set up the run script for Nenzfr.
    #paramDict = {'jobName':quote(opt.jobName), 'gasName':quote(opt.gasName),
    #             'T1':opt.T1, 'p1':opt.p1, 'Vs':opt.Vs, 'pe':opt.pe,
    #             'contourFileName':quote(opt.contourFileName),
    #             'gridFileName':quote(opt.gridFileName), 
    #             'chemModel':quote(opt.chemModel),
    #             'areaRatio':opt.areaRatio,
    #             'nni':opt.nni, 'nnj':opt.nnj, 'nbi':opt.nbi, 'bx':opt.bx, 'by':opt.by,
    #             'max_time':opt.max_time, 'max_step':opt.max_step}
    #prepare_input_script(paramDict, opt.jobName)


    
    # Run Eilmer3
    #run_command(E3BIN+('/e3prep.py --job=%s --do-svg' % (opt.jobName,)))
    #run_command(E3BIN+('/e3shared.exe --job=%s --run' % (opt.jobName,)))
    
     
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
