#!/usr/bin/env python
# run_adaptive_simulation.py
#
# Top-level script to coordinate the running of the 
# solution-adaptive simulation in stages.
#
# PJ, 22-Feb-2010, 19-Apr-2010

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
sys.path.append(os.path.expandvars("$HOME/e3bin"))

#---------------------------------------------------------------
def prepare_input_script(substituteDict, jobName, stage):
    """
    Prepare the actual input file from a template.
    """
    stageName = jobName + str(stage)
    templateFileName = jobName + ".template"
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
    installDir = "/home/paul/e3bin"
    run_command("%s/e3prep.py --job=%s --do-svg" % (installDir, stageName,))
    run_command("%s/e3shared.exe --job=%s --run" % (installDir, stageName,))
    return

#---------------------------------------------------------------
# main...

jobName = 'nasaRIT'

lowResStatorMesh = [{'i':12,'j':10,'k':10}, {'i':12,'j':10,'k':10}, {'i':4,'j':10,'k':10}]
lowResRotorMesh = [{'i':4,'j':10,'k':10}, {'i':30,'j':10,'k':10}, {'i':20,'j':10,'k':10}]
lowResBeta = 1.08
medResStatorMesh = [{'i':17,'j':15,'k':15}, {'i':17,'j':15,'k':15}, {'i':5,'j':15,'k':15}]
medResRotorMesh = [{'i':5,'j':15,'k':15}, {'i':49,'j':15,'k':15}, {'i':30,'j':15,'k':15}]
medResBeta = 1.07
highResStatorMesh = [{'i':17,'j':30,'k':30}, {'i':17,'j':30,'k':30}, {'i':5,'j':30,'k':30}]
highResRotorMesh = [{'i':5,'j':30,'k':30}, {'i':49,'j':30,'k':30}, {'i':30,'j':30,'k':30}]
highResBeta = 1.07

stage = 0
paramDict = {'jobName': jobName, 
             'stage': stage,
             'stator_nijk': highResStatorMesh,
             'rotor_nijk': highResRotorMesh,
             'diffuse_beta': highResBeta, 
             'cfl':0.5, 
             'sameGrid':1, 
             'dt_init':1.0e-7, 
             'max_step':100}
run_stage(paramDict, jobName, stage)


print "Done at top-level."

