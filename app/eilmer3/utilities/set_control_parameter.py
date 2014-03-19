#!/usr/bin/env python

"""
A Simple python script to set a parameter in an eilmer3 .control file.
"""

import sys

if __name__=="__main__":
    if len(sys.argv)!=4:
        print "Usage: set_control_parameter.py <controlFile> <parameterName> <value>"
        sys.exit()
    controlFile = sys.argv[1]
    parameterName = sys.argv[2]
    value = sys.argv[3]

    ifile = open(controlFile,"r")
    lines = ifile.readlines()
    ifile.close()

    newLines = ""

    found = False
    for i,line in enumerate(lines):
        tks = line.split("=")
        print tks
        if tks[0].replace(" ","")==parameterName:
            newLines += "%s = %s\n" % ( parameterName, value )
            found = True
        else:
            newLines.append(line)

    if not found:
        print "parameterName: %s was not found in file: %s" % ( parameterName, controlFile )
        sys.exit()

    # os.system("mv %s %s.org" % ( controlFile, controlFile ) )
    
    ofile = open(controlFile,"w")
    ofile.write(newLines)
    ofile.close()

    print "Done."
