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

VERSION_STRING = "16-April-2012"

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
import optparse
#import numpy
from numpy import array, mean, logical_and, zeros
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#---------------------------------------------------------------

def prepare_run_script(substituteDict, jobName, Cluster):
    """
    Prepare the actual run file for Nenzfr from a template.
    """
    templateFileName = "run_template_"+Cluster +".sh"
    # Check that the templateFileName is in the current directory,
    # if it isn't copy it from "E3BIN/nenzfr_data_files/"
    if not os.path.exists(templateFileName):
        command_text = 'cp '+E3BIN+'/nenzfr_data_files/'+templateFileName+' ./'
        #print command_text
        run_command('cp '+E3BIN+'/nenzfr_data_files/'+templateFileName+' ./')
    scriptFileName = "run_" + jobName + ".sh"
    fp = open(templateFileName, 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.safe_substitute(substituteDict)
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
        #print "k= %d; tstart= %f; pe[k]= %f" % (k, tstart+k*dt, pe[k])
    
    return pe


def main():
    """
    Examine the command-line options to decide the what to do
    and then coordinate the calculations done by Nenzfr.
    """
    op = optparse.OptionParser(version=VERSION_STRING)

    op.add_option('--runCMD', dest='runCMD', default='./',
                  choices=['./', 'qsub '],
                  help=("command used to execute the run-script file "
                        "[default: %default]"))
    op.add_option('--Cluster', dest='Cluster', default='Mango',
                  choices =['Mango', 'Barrine'],
                  help=("specify on which cluster the computations are to be ran. "
                        "This is used to define which run template script will "
                        "be used. Options: "
                        "Mango; Barrine [default: %default]"))

    op.add_option('--tstart', dest='tstart', type='float', default='1.5e-3',
                  help=("time at which the first slice of the input pressure "
                        "trace is to begin [default: %default]"))
    op.add_option('--nsteps', dest='nsteps', type='int', default=5,
                  help=("number of slices to use [default: %default]"))
    op.add_option('--dt', dest='dt', type='float', default=0.5e-3,
                  help=("width of each averaging slice [default: %default]"))

    op.add_option('--pefile', dest='peFileName', default=None,
                  help="file specifying the transient equilibrium pressure (in Pa) "
                       "[default: %default]")
    op.add_option('--job', dest='jobName', default='nozzle',
                  help="base name for Eilmer3 files [default: %default]")

    # Multitude of options required by nenzfr.
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
    op.add_option('--chem', dest='chemModel', default='eq',
                  choices=['eq', 'neq', 'frz', 'frz2'],
                  help=("chemistry model: " "eq=equilibrium; "
                        "neq=non-equilibrium; " "frz=frozen "
                        "[default: %default]"))
    op.add_option('--area', dest='areaRatio', default=1581.165,
                  help=("nozzle area ratio. only used for estcj calc. "
                        "use when --cfile(--gfile) are "
                        "specified. [default: %default]"))
    op.add_option('--cfile', dest='contourFileName',
                  default='contour-t4-m10.data',
                  help="file containing nozzle contour [default: %default]")
    op.add_option('--gfile', dest='gridFileName', default='None',
                  help="file containing nozzle grid. "
                  "overrides --cfile if both are given "
                  "[default: %default]")
    op.add_option('--exitfile', dest='exitSliceFileName',
                  default='nozzle-exit.data',
                  help="file for holding the nozzle-exit data [default: %default]")
    op.add_option('--block-marching', dest='blockMarching', action='store_true',
                  default=False, help="run nenzfr in block-marching mode")
    # The following defaults suit a Mach 10 Nozzle calculation.
    op.add_option('--nni', dest='nni', type='int', default=1800,
                  help=("number of axial cells"))
    op.add_option('--nnj', dest='nnj', type='int', default=100,
                  help=("number of radial cells"))
    op.add_option('--nbi', dest='nbi', type='int', default=180,
                  help=("number of axial blocks for the divergence section (nozzle_blk)"))
    op.add_option('--nbj', dest='nbj', type='int', default=1,
                  help=("number of radial blocks"))
    op.add_option('--bx', dest='bx', type='float', default=1.05,
                  help=("clustering in the axial direction"))
    op.add_option('--by', dest='by', type='float', default=1.002,
                  help=("clustering in the radial direction"))
    op.add_option('--max-time', dest='max_time', type='float', default=6.0e-3,
                  help=("overall simulation time for nozzle flow"))
    op.add_option('--max-step', dest='max_step', type='int', default=80000,
                  help=("maximum simulation steps allowed"))

    op.add_option('--Twall', dest='Tw', type='float', default=300.0,
                  help=("Nozzle wall temperature, in K "
                        "[default: %default]"))
    op.add_option('--BLTrans', dest='BLTrans', default="x_c[-1]*1.1",
                  help=("Transition location for the Boundary layer. Used "
                        "to define the turbulent portion of the nozzle. "
                        "[default: >nozzle length i.e. laminar nozzle]"))
    op.add_option('--TurbVisRatio', dest='TurbVisRatio', type='float',
                  default=100.0, help=("Turbulent to Laminar Viscosity Ratio "
                  "[default: %default]"))
    op.add_option('--TurbIntensity', dest='TurbInten', type='float',
                  default=0.05, help=("Turbulence intensity at the throat "
                  "[default: %default]"))
    op.add_option('--CoreRadiusFraction', dest="coreRfraction", type='float',
                  default=2.0/3.0, help=("Radius of core flow as a fraction of "
                  "the nozzle exit radius [default: %default]"))
    opt, args = op.parse_args()
    
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
    if opt.blockMarching is True:
        opt.blockMarching = "--block-marching"
    else:
        opt.blockMarching = ""
    if bad_input:
        return -2
    
    # Calculate the averaged equilibrium pressure for each
    # time interval...
    pe = calculate_pe_values(opt.peFileName, opt.tstart, opt.nsteps, opt.dt)
    
    # Set up the run script for Nenzfr.
    paramDict = {'jobName':quote(opt.jobName), 'gasName':quote(opt.gasName),
                 'T1':opt.T1, 'p1':opt.p1, 'Vs':opt.Vs,
                 'chemModel':quote(opt.chemModel),
                 'contourFileName':quote(opt.contourFileName),
                 'gridFileName':quote(opt.gridFileName),
                 'exitSliceFileName':quote(opt.exitSliceFileName),
                 'areaRatio':opt.areaRatio, 'blockMarching':opt.blockMarching,
                 'nni':opt.nni, 'nnj':opt.nnj, 'nbi':opt.nbi, 'nbj':opt.nbj,
                 'bx':opt.bx, 'by':opt.by,
                 'max_time':opt.max_time, 'max_step':opt.max_step,
                 'Tw':opt.Tw, 'TurbVisRatio':opt.TurbVisRatio,
                 'TurbInten':opt.TurbInten, 'BLTrans':opt.BLTrans,
                 'CoreRadiusFraction':opt.coreRfraction}
    
    # As building an equilibrium gas LUT is so time consuming, we do it here
    # and then copy the resulting LUT into each case sub-directory. The following
    # lines are copied almost verbatim from "nenzfr.py"
    if opt.chemModel in ['eq']:
        if opt.gasName in ['n2']:
            eqGasModelFile = 'cea-lut-'+upper(opt.gasName)+'.lua.gz'
        else:
            eqGasModelFile = 'cea-lut-'+opt.gasName+'.lua.gz'
        if not os.path.exists(eqGasModelFile):
            run_command('build-cea-lut.py --gas='+opt.gasName)
        paramDict['gmodelFile'] = eqGasModelFile

    
    for k in range(len(pe)):
        # Create sub-directory for the current case.
        caseString = 'case'+"{0:03}".format(k)
        
        print 60*"-"
        print caseString
        print "tstart= %f; pe[k]= %f" % (opt.tstart+k*opt.dt, pe[k])
        
        run_command('mkdir ./'+caseString)
        
        # Set up the run script for Nenzfr
        #paramDict = {'jobName':quote(opt.jobName), 'gasName':quote(opt.gasName),
        #             'T1':opt.T1, 'p1':opt.p1, 'Vs':opt.Vs, 'pe':pe[k],
        #             'chemModel':quote(opt.chemModel),
        #             'contourFileName':quote(opt.contourFileName),
        #             'gridFileName':quote(opt.gridFileName),
        #             'exitSliceFileName':quote(opt.exitSliceFileName),
        #             'areaRatio':opt.areaRatio, 'blockMarching':opt.blockMarching,
        #             'nni':opt.nni, 'nnj':opt.nnj, 'nbi':opt.nbi, 'nbj':opt.nbj,
        #             'bx':opt.bx, 'by':opt.by,
        #             'max_time':opt.max_time, 'max_step':opt.max_step,
        #             'Tw':opt.Tw, 'TurbVisRatio':opt.TurbVisRatio,
        #             'TurbInten':opt.TurbInten, 'BLTrans':opt.BLTrans,
        #             'CoreRadiusFraction':opt.coreRfraction}
        paramDict['pe'] = pe[k]
        scriptFileName = prepare_run_script(paramDict, opt.jobName+'_'+caseString, opt.Cluster)
         
        # Move the run script to its sub-directory
        command_text = 'mv '+scriptFileName+' ./'+caseString+'/'+scriptFileName
        run_command(command_text)
        
        # If required, copy the nozzle.timing file to the sub-directory
        if opt.blockMarching is True:
            if os.path.exists(nozzle.timing):
                command_text = 'cp nozzle.timing ./'+caseString+'/'
                run_command(command_text)
        
        # If required, copy the nozzle.timing file to the sub-directory
        if paramDict['blockMarching'] in ["--block-marching",]:
            if os.path.exists('nozzle.timing'):
                command_text = 'cp nozzle.timing ./'+caseString+'/'
                run_command(command_text)

        # If require, copy the equilibrium gas LUT to the sub-directory
        if paramDict['chemModel'] in ['"eq"',]:
            command_text = 'cp '+paramDict['gmodelFile']+' ./'+caseString+'/'
            run_command(command_text)
     
        # Change into the sub-directory, ensure the run script is exectuable and
        # then run it
        os.chdir(caseString)
        run_command('chmod u+x '+scriptFileName)
        print ""
        print opt.runCMD+scriptFileName
        print ""
        # I am not sure how to replace the next line with the run_command function
        os.system(opt.runCMD+scriptFileName)
        os.chdir('../')
    
    return 0

#---------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr Quasi Transient:\n Calculate Shock Tunnel Test Flow Conditions for a varying Supply Pressure"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
