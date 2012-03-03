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

VERSION_STRING = "02-Mar-2012"

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

    op.add_option('--runCMD', dest='runCMD', default='qsub ',
                  help=("command used to execute the run-script file "
                        "[default: %default]"))
    op.add_option('--Cluster', dest='Cluster', default='Barrine',
                  choices =['Mango', 'Barrine'],
                  help=("specify on which cluster the computations are to be ran. "
                        "This is also used to define which run template script will "
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
    opt, args = op.parse_args()
    
    # Go ahead with a new calculation.
    # First, make sure that we have the needed parameters.
    bad_input = False
    
    if opt.peFileName is None:
        print "Need to supply a string value for pefile."
        bad_input = True
    if bad_input:
        return -2
    
    # Calculate the averaged equilibrium pressure for each
    # time interval...
    pe = calculate_pe_values(opt.peFileName, opt.tstart, opt.nsteps, opt.dt)
    
    # Set up the run script for Nenzfr.
    paramDict = {'jobName':quote(opt.jobName)}
    for k in range(len(pe)):
        # Create sub-directory for the current case.
        caseString = 'case'+"{0:02}".format(k)
        
        print 60*"-"
        print caseString
        print "tstart= %f; pe[k]= %f" % (opt.tstart+k*opt.dt, pe[k])
        
        run_command('mkdir ./'+caseString)
        
        # Set up the run script for Nenzfr
        paramDict['pe'] = pe[k]
        scriptFileName = prepare_run_script(paramDict, opt.jobName+'_'+caseString, opt.Cluster)
         
        # Move the run script to its sub-directory
        command_text = 'mv '+scriptFileName+' ./'+caseString+'/'+scriptFileName
        run_command(command_text)
        
        # Change into the sub-directory, ensure the run script is exectuable and
        # then run it
        os.chdir(caseString)
        run_command('chmod u+x '+scriptFileName)
        print ""
        print opt.runCMD+scriptFileName
        print ""
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
