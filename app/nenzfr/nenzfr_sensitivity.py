#!/usr/bin/env python
# nenzfr_sensitivity.py
#
# This script coordinates the running of a series of nenzfr 
# calculations... 
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


def configure_input_dictionary(perturbedDict, defaultPerturbations, gradient, perturb): 
    if perturbedDict["BLTrans"] == "x_c[-1]*1.1":
        del perturbedDict["BLTrans"]
    for i in range(len(perturbedDict)):
            var = perturbedDict.keys()[i]
            temp =  perturbedDict.values()[i]
            temp = temp.strip(']').strip('[').split(',')
            temp = [float(temp[j]) for j in range(len(temp))]
            perturbedDict[var] = set_perturbed_values(var, temp, defaultPerturbations, perturb, gradient) 
    return perturbedDict

def set_perturbed_values(var, temp, defaultPerturbations, perturb, gradient):
    if gradient == "linear":
        if len(temp) == 1:
            # have just the value
            print "Using default relative perturbation for "+var
            temp = [temp[0], temp[0]*(1+defaultPerturbations[var]), temp[0]*(1-defaultPerturbations[var])]
        elif len(temp) == 2:
            # have either [value, delta] or [value, percent_delta]
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]), temp[0]*(1.0-temp[1])]
            else:
                temp = [temp[0], temp[0]+temp[1], temp[0]-temp[1]]
        elif len(temp) == 3:
            # have either [value, +delta, -delta] or [value, +percent_delta, -percent_delta]
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]), temp[0]*(1.0+temp[2])]
            else:
                temp = [temp[0], temp[0]+temp[1], temp[0]+temp[2]]
        else:
            print "Too many values given for "+var
            return -2
    elif gradient == "quad":
        if len(temp) == 1:
            print "Using default relative perturbation for "+var
            temp = [temp[0], temp[0]*(1+defaultPerturbations[var]), temp[0]*(1-defaultPerturbations[var]),\
                    temp[0]*(1+defaultPerturbations[var]*2.0), temp[0]*(1-defaultPerturbations[var]*2.0)]
        elif len(temp) == 2:
            # have either [value, delta] or [value, +percent_delta]
            # we assume that "2"*delta == 2*delta
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]), temp[0]*(1.0-temp[1]),\
                        temp[0]*(1.0+temp[1]*2.0), temp[0]*(1.0-temp[1]*2.0)]
            else:
                temp = [temp[0], temp[0]+temp[1], temp[0]-temp[1],\
                        temp[0]+2.0*temp[1], temp[0]-2.0*temp[1]]
        elif len(temp) == 3:
            # have either [value, delta, "2"*delta] or [value, percent_delta, "2"*percent_delta]
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]), temp[0]*(1.0-temp[1]),\
                        temp[0]*(1.0+temp[2]), temp[0]*(1.0-temp[2])]
            else:
                temp = [temp[0], temp[0]+temp[1], temp[0]-temp[1],\
                        temp[0]+temp[2], temp[0]-temp[2]]
        elif len(temp) == 4:
            # have either [value, +delta, -delta, "2"*delta] or
            #             [value, +percent_delta, -percent_delta, "2"*percent_delta]
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]), temp[0]*(1.0-temp[2]),\
                        temp[0]*(1.0+temp[3]), temp[0]*(1.0-temp[3])]
            else:
                temp = [temp[0], temp[0]+temp[1], temp[0]-temp[2],\
                        temp[0]+temp[3], temp[0]-temp[3]]
        elif len(temp) == 5:
            # have either [value, +delta, -delta, +"2"*delta, -"2"*delta] or
            #             [value +percent_delta, -percent_delta, +"2"*percent_delta, -"2"*percent_delta]
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]), temp[0]*(1.0-temp[2]),\
                        temp[0]*(1.0+temp[3]), temp[0]*(1.0-temp[4])]
            else:
                temp = [temp[0], temp[0]+temp[1], temp[0]-temp[2],\
                        temp[0]+temp[3], temp[0]-temp[4]]
        else:
            print "Too many values given for "+var
            return -2
    return temp


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
                        "This is used to define which run template script will "
                        "be used. Options: "
                        "Mango; Barrine [default: %default]"))
    
    op.add_option('--perturb', dest='perturb', default='rel',
                  help=("specify whether the given perturbations are absolute "
                        "or relative values. Options: rel, abs [default: %default]"))

    op.add_option('--gradient', dest='gradient', default='linear',
                  help=("specify whether the gradient is to be calculated "
                        "using a linear or quadratic equation. ONLY available "
                        "with the --create-LUT option "
                        "Options: linear, quadratic [default: %default]"))

    op.add_option('--job', dest='jobName', default='nozzle',
                  help="base name for Eilmer3 files [default: %default]")

    # Multitude of options required by nenzfr.
    op.add_option('--gas', dest='gasName', default='air',
                  choices=['air', 'air5species', 'n2', 'co2', 'h2ne'],
                  help=("name of gas model: "
                        "air; " "air5species; " "n2; " "co2; " "h2ne"))
    #op.add_option('--p1', dest='p1', type='float', default=None,
    #              help=("shock tube fill pressure or static pressure, in Pa"))
    
    op.add_option('--p1', dest='p1', default=None, 
                  help=("shock tube fill pressure (in Pa) and its perturbation"))
    
    op.add_option('--T1', dest='T1', default=None,
                  help=("shock tube fill temperature, in degrees K"))
    op.add_option('--Vs', dest='Vs', default=None,
                  help=("incident shock speed, in m/s"))
    op.add_option('--pe', dest='pe', default=None,
                  help=("equilibrium pressure (after shock reflection), in Pa"))
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
    
    op.add_option('--create-LUT', dest='createLUT', action='store_true',
                  default=False, 
                  help="create a LUT by perturbing only p1, T1, Vs and pe")

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

    op.add_option('--Twall', dest='Tw', default='300.0',
                  help=("Nozzle wall temperature, in K "
                        "[default: %default]"))
    op.add_option('--BLTrans', dest='BLTrans', default="x_c[-1]*1.1",
                  help=("Transition location for the Boundary layer. Used "
                        "to define the turbulent portion of the nozzle. "
                        "[default: >nozzle length i.e. laminar nozzle]"))
    op.add_option('--TurbVisRatio', dest='TurbVisRatio', default='100.0', 
                  help=("Turbulent to Laminar Viscosity Ratio "
                  "[default: %default]"))
    op.add_option('--TurbIntensity', dest='TurbInten', default='0.05', 
                  help=("Turbulence intensity at the throat "
                  "[default: %default]"))
    op.add_option('--CoreRadiusFraction', dest="coreRfraction", default='0.6666666667', 
                  help=("Radius of core flow as a fraction of "
                  "the nozzle exit radius [default: %default]"))
    opt, args = op.parse_args()
    
    
    # Set all the default relative perturbation values
    defaultPerturbations = {'p1':0.05, 'T1':0.05, 'Vs':0.05, 'pe':0.05, 
                           'Tw':0.05, 'BLTrans':0.05, 'TurbVisRatio':0.05,
                           'TurbInten':0.05, 'coreRfraction':0.05}

    # Go ahead with a new calculation.
    # First, make sure that we have the needed parameters.
    bad_input = False

    #print opt.p1[2]
    #test = opt.p1.strip('[').strip(']').split(',')
    #print test
    #for i in range(len(test)):
    #    test[i] = float(test[i])
    #print test
    #test2[i] = float(test2[i])
    #p1 = {}
    #for i in range(len(items)):
    #    p1[variable_list[i]].append(float(items[i]))
    
    #print opt
    #print opt['p1']
    #print opt.values()
    #print opt.values().strip('[').strip(']').split(',')
    
    if opt.p1 is None:
        print "Need to supply a value for p1."
        bad_input = True
    if opt.T1 is None:
        print "Need to supply a value for T1."
        bad_input = True
    if opt.Vs is None:
        print "Need to supply a value for Vs."
        bad_input = True    
    if opt.pe is None:
        print "Need to supply a value for pe."
        bad_input = True
    if opt.blockMarching is True:
        opt.blockMarching = "--block-marching"
    else:
        opt.blockMarching = ""
    if opt.createLUT is False:
        if opt.gradient is "quadratic":
            opt.gradient = 'linear'
            print "Sensitivity will be calculated using a linear slope."
    if bad_input:
        return -2
    
    #print opt
    #print opt.strip('[').strip(']').split(',')
    if opt.createLUT:
        #perturbedVariables = ['p1','T1','Vs','pe']
        perturbedDict = {'p1':opt.p1, 'T1':opt.T1, 'Vs':opt.Vs, 'pe':opt.pe}
        perturbedDict = configure_input_dictionary(perturbedDict, defaultPerturbations, 
                                                opt.gradient, opt.perturb)
    else:
        #perturbedVariables = ['p1','T1','Vs','pe','Tw','BLTrans','TurbVisRatio',
        #                      'TurbInten','coreRfraction',]
        perturbedDict = {'p1':opt.p1, 'T1':opt.T1, 'Vs':opt.Vs, 'pe':opt.pe,
                        'Tw':opt.Tw, 'BLTrans':opt.BLTrans, 'TurbVisRatio':opt.TurbVisRatio,
                        'TurbInten':opt.TurbInten, 'coreRfraction':opt.coreRfraction}
        perturbedDict = configure_input_dictionary(perturbedDict, defaultPerturbations,
                                                opt.gradient, opt.perturb)

    #inputDict = {'p1':[float(opt.p1.strip('[').strip(']').split(',')[i]) for i in range(len(a))]
    
    print perturbedDict
    
    # Set up the run script for Nenzfr.
    paramDict = {'jobName':quote(opt.jobName)}
    #for k in range(len(pe)):
    #    # Create sub-directory for the current case.
    #    caseString = 'case'+"{0:03}".format(k)
    #    
    #    print 60*"-"
    #    print caseString
    #    print "tstart= %f; pe[k]= %f" % (opt.tstart+k*opt.dt, pe[k])
    #    
    #    run_command('mkdir ./'+caseString)
    #    
    #    # Set up the run script for Nenzfr
    #    paramDict = {'jobName':quote(opt.jobName), 'gasName':quote(opt.gasName),
    #                 'T1':opt.T1, 'p1':opt.p1, 'Vs':opt.Vs, 'pe':pe[k],
    #                 'chemModel':quote(opt.chemModel),
    #                 'contourFileName':quote(opt.contourFileName),
    #                 'gridFileName':quote(opt.gridFileName),
    #                 'exitSliceFileName':quote(opt.exitSliceFileName),
    #                 'areaRatio':opt.areaRatio, 'blockMarching':opt.blockMarching,
    #                 'nni':opt.nni, 'nnj':opt.nnj, 'nbi':opt.nbi, 'nbj':opt.nbj,
    #                 'bx':opt.bx, 'by':opt.by,
    #                 'max_time':opt.max_time, 'max_step':opt.max_step,
    #                 'Tw':opt.Tw, 'TurbVisRatio':opt.TurbVisRatio,
    #                 'TurbInten':opt.TurbInten, 'BLTrans':opt.BLTrans,
    #                 'CoreRadiusFraction':opt.coreRfraction}
    #    scriptFileName = prepare_run_script(paramDict, opt.jobName+'_'+caseString, opt.Cluster)
    #     
    #    # Move the run script to its sub-directory
    #    command_text = 'mv '+scriptFileName+' ./'+caseString+'/'+scriptFileName
    #    run_command(command_text)
    #    
    #    # Change into the sub-directory, ensure the run script is exectuable and
    #    # then run it
    #    os.chdir(caseString)
    #    run_command('chmod u+x '+scriptFileName)
    #    print ""
    #    print opt.runCMD+scriptFileName
    #    print ""
    #    os.system(opt.runCMD+scriptFileName) # I am not sure how to replace this with the run_command function
    #    os.chdir('../')
     
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
