#! /usr/bin/python

import sys
import os
from getopt import getopt, GetoptError
from gaspy import *

longOptions = ["help", "model=", "species=", "output="]

def printUsage():
    print ""
    print "Usage: gasmodel.py [--help] [--model=<modelName>] [--species=<speciesList|none>]"
    print "                   [--output=<luaFile|gas-model.lua>]\n"
    print "e.g. gasmodel.py --model='thermally perfect gas' --species='N2 N' --output='nitrogen.lua'\n"
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

    create_gas_file( model, species, output )
    print "Gas model created in file:", output
   
    print "done."
    
if __name__ == '__main__':
    main()
