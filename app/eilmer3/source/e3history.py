#! /usr/bin/env python
"""
e3history.py -- Pick up history files and extract data for a specified cell.

e3history.py is one of the eilmer3 postprocessing programs.  
Try invoking it with the --help option to see more information.

.. Author: P.Jacobs

.. Versions:
   19-June-2008
   01-Sep-2008 added dictionary of data values and Pitot pressure calc.
"""

# ----------------------------------------------------------------------
#
import sys
import os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from getopt import getopt
from libprep3 import create_gas_model

shortOptions = ""
longOptions = ["help", "ijk=", "cell=", "output=", "add-pitot-p", "add-total-p",
               "add-total-enthalpy", "add-mach", "add-molef", "gmodel-file="]

def printUsage():
    print ""
    print "Usage: e3history.py" + \
          " [--help]" + \
          " [--ijk=<i,j,k>]" + \
          " [--cell=<index>]" + \
          " [--output=<filename>]" + \
          " [--add-pitot-p]" + \
          " [--add-total-p]" + \
          " [--add-total-enthalpy]" + \
          " [--add-mach]" + \
          " [--add-molef --gmodel-file=\"gas-model.lua\"]" + \
          " <history-file>"
    print ""
    print "Notes:"
    print "  If you just specify the history-file name, you will get a summary of its content."
    print "  The cell order in the history file starts at 1."
    print "  Cell indexing within a block grid starts at 0."
    print "  History files are stored in the hist/ folder."
    print ""
    print "Examples:"
    print "  $ ./e3history.py --ijk=\"358,4,0\" --add-pitot-p hist/m4cone.hist.b0000"
    print "  $ ./e3history.py --cell=3 --output=test2.txt  hist/m4cone.hist.b0000"
    return

def get_variable_list(line):
    tokens = line.split()
    del tokens[0]
    var_list = []
    for token in tokens:
        var_list.append(token.strip('"')) # just keep the name
    return var_list

#----------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin e3history.py..."
    #
    if sys.argv[-1] == "--help":
        printUsage()
        print "Done."
        sys.exit(0)
    userOptions = getopt(sys.argv[1:-1], shortOptions, longOptions)
    if len(sys.argv) > 1: 
        historyFile = sys.argv[-1]
    else:
        historyFile = ""
    uoDict = dict(userOptions[0])
    if (len(uoDict)==0 and len(historyFile)==0) or uoDict.has_key("--help"):
        printUsage()
        print "Done."
        sys.exit(0)
    if not os.path.exists(historyFile):
        print "Error: could not find history file:", historyFile
        print "Done."
        sys.exit(-1)
    if uoDict.has_key("--ijk"):
        ijkString = uoDict.get("--ijk", "0,0,0")
        tokens = ijkString.split(',')
        iSelect = int(tokens[0])
        jSelect = int(tokens[1])
        kSelect = int(tokens[2])
        print "Pull out cell with i=", iSelect, "j=", jSelect, "k=", kSelect, "indices."
    if uoDict.has_key("--cell"):
        iCell = int(uoDict.get("--cell", 0))
        print "Pull out cell number", iCell
    # Always try to attach a gas model, possibly using the command-line argument
    # or by falling back to the expected default file name.
    gmodelFileName = uoDict.get("--gmodel-file", "gas-model.lua")
    print "Attempt to create gas model from file:", gmodelFileName
    if os.path.exists(gmodelFileName):
        gmodel = create_gas_model(gmodelFileName)
        nsp = gmodel.get_number_of_species()
        speciesList = [gmodel.species_name(isp) for isp in range(nsp)]
        print "speciesList=", speciesList
    else:
        print "Failed to create gas model."
        gmodel = None
        nsp = 1
        speciesList = []
    # Sift through the command-line dictionary
    # to determine which extra variables should be added.
    add_pitot_p = uoDict.has_key("--add-pitot-p")
    add_total_p = uoDict.has_key("--add-total-p")
    add_total_enthalpy = uoDict.has_key("--add-total-enthalpy")
    add_mach = uoDict.has_key("--add-mach")
    add_molef = uoDict.has_key("--add-molef") and (gmodel != None)
    #
    # If we reach this point, go ahead and try to analyse its content and,
    # maybe, output a selection.
    #
    print "History File:", historyFile
    outputFile = uoDict.get("--output", "output.txt")
    fp_in = open(historyFile, "r")
    fp_out = open(outputFile, "w")
    # Expect the first line to list the variables on each data line.
    line = fp_in.readline()
    var_list = get_variable_list(line)
    print "Have read", len(var_list), "variables, specifically", var_list

    fp_out.write("#")  # first line gives variable names as a comment to GNUplot
    var_count = 0
    for var_name in var_list:
        var_count += 1
        fp_out.write(" %d:%s" % (var_count, var_name) )
    if add_pitot_p:
        var_count += 1
        fp_out.write(" %d:pitot_p" % var_count)
    if add_total_p:
        var_count += 1
        fp_out.write(" %d:total_p" % var_count)
    if add_total_enthalpy:
        var_count += 1
        fp_out.write(" %d:total_h" % var_count)
    if add_mach:
        var_count += 1
        fp_out.write(" %d:mach" % var_count)
    if add_molef:
        for isp in range(nsp):
            specname = gmodel.species_name(isp).replace(' ', '-')
            varName = "molef[%d]-%s" % (isp, specname)
            var_count += 1
            fp_out.write(" %d:%s" % (var_count, varName))
    fp_out.write("\n")

    # At this point, we expect the line to contain significant data.
    firstDataLine = True
    while 1:
        line = fp_in.readline()
        if len(line) == 0: break
        if line[0] == '#': continue  # ignore comments in the file
        tokens = line.split()
        # Build a dictionary or the variables and their associated values
        data = {}
        for i in range(len(var_list)):
            data[var_list[i]] = tokens[i]
        if firstDataLine:
            oldTimeValue = float(data['time'])
            cellCount = 1
            timeCount = 1
            firstDataLine = False
        else:
            newTimeValue = float(data['time'])
            if abs(newTimeValue - oldTimeValue) < 1.0e-10:
                cellCount += 1
            else:
                oldTimeValue = newTimeValue
                timeCount += 1
                cellCount = 1
        if (uoDict.has_key("--cell") and cellCount==iCell) or \
                (uoDict.has_key("--ijk") and int(data['i'])==iSelect and \
                     int(data['j'])==jSelect and int(data['k'])==kSelect):
            # Write out the cell's original data
            for i in range(len(var_list)):
                fp_out.write(" %s" % tokens[i])
            # and, maybe, add Pitot pressure, total pressure and Mach number.
            ux = float(data['vel.x'])
            uy = float(data['vel.y'])
            uz = float(data['vel.z'])
            a = float(data['a'])
            p = float(data['p'])
            rho = float(data['rho'])
            e0 = float(data['e[0]'])
            u = math.sqrt(ux*ux + uy*uy + uz*uz)
            M = u/a
            tke = float(data['tke'])
            g = a*a*rho/p # approximation for gamma
            if add_pitot_p:
                # Rayleigh Pitot formula
                if M > 1.0:
                    # Go through the shock and isentropic compression.
                    t1 = (g+1)*M*M/2;
                    t2 = (g+1)/(2*g*M*M - (g-1));
                    pitot_p = p * math.pow(t1,(g/(g-1))) * math.pow(t2,(1/(g-1)));
                else:
                    # Isentropic compression only.
                    t1 = 1 + 0.5*(g-1)*M*M;
                    pitot_p = p * math.pow(t1,(g/(g-1)));
                fp_out.write(" %e" % pitot_p)
            if add_total_p:
                # Isentropic compression only.
                t1 = 1 + 0.5*(g-1)*M*M;
                total_p = p * math.pow(t1,(g/(g-1)));
                fp_out.write(" %e" % total_p)
            if add_total_enthalpy:
                # Sum up the bits of energy,
                # forgetting the multiple energy modes, for the moment.
                total_h = p/rho + e0 + 0.5*u*u + tke;
                fp_out.write(" %e" % total_h)
            if add_mach:
                fp_out.write(" %e" % M)
            if add_molef:
                # Accumulate the mass fractions and then compute the mole fractions
                # in the context of the gas model.
                massf = []
                for isp in range(nsp):
                    specname = gmodel.species_name(isp).replace(' ', '-')
                    massf.append(float(data["massf[%d]-%s" % (isp, specname)]))
                molef = gmodel.to_molef(massf)
                for isp in range(nsp):
                    fp_out.write(" %e" % molef[isp])
            # We have finished writing values; terminate the output line.
            fp_out.write("\n") 

    print "cellCount=", cellCount, "timeCount=", timeCount
    fp_in.close()
    fp_out.close()

    print "Done."

