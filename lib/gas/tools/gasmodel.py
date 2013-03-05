#! /usr/bin/python
# File: gasmodel.py
# Author: Daniel Potter
# Notes: A python tool to create a gas model file from the command line

import sys
import os
from getopt import getopt, GetoptError
from gaspy import *

longOptions = ["help", "model=", "species=", "output=", "lut-file="]

def printUsage():
    print ""
    print "Usage: gasmodel.py [--help] [--model=<modelName>] [--species=<speciesList|none>]"
    print "                   [--lut-file=<LUTFileName>]"
    print "                   [--output=<luaFile|gas-model.lua>]"
    print ""
    print "Examples:"
    print "$ gasmodel.py --model='thermally perfect gas' --species='N2 N' --output='nitrogen.lua'"
    print ""
    print "$ gasmodel.py --model='ideal gas' --species='Ar He' --lut-file='cea-lut-custom.lua.gz' \\"
    print "              --output='LUT-plus-Ar-He.lua'"
    print ""
    return

def main():
    #
    try:
        userOptions = getopt(sys.argv[1:], [], longOptions)
    except GetoptError, e:
        print "One (or more) of your command-line options was no good."
        print "    ", e
        printUsage()
        sys.exit(1)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        printUsage()
        sys.exit(0)
    #
    model = uoDict.get("--model", "none")
    species = uoDict.get("--species", "").split()
    print "species list:", species
    output = uoDict.get("--output","gas-model.lua")
    lut_file = uoDict.get("--lut-file", None)
    if lut_file:
        print "Will prepend LUT gas species from", lut_file
    #
    create_gas_file(model, species, output, lut_file)
    print "Gas model created in file:", output
    #
    print "done."
    return
    
if __name__ == '__main__':
    main()
