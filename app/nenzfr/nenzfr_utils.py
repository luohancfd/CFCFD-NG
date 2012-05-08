"""
nenzfr_utils.py -- Small utility functions needed by the main program.
"""

import sys, os
import shlex, subprocess, string
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
    if (type(cmdText) is list):
        args = cmdText
    else:
        args = shlex.split(cmdText)
    print "About to run cmd:", string.join(args)
    # p = subprocess.Popen(args)
    # wait until the subprocess is finished
    # stdoutData, stderrData = p.communicate()
    return subprocess.check_call(args)

def quote(str):
    """
    Put quotes around a string.
    """
    return '"' + str + '"'

