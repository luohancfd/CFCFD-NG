#! /usr/bin/python
"""
A python tool to create a gas model file from the command line.

Try the command::

  gasmodel.py --help

..File: gasmodel.py

..Author: Daniel Potter (with docs by PJ)
"""

import sys
import os
from getopt import getopt, GetoptError
from gaspy import *

longOptions = ["help", "model=", "species=", "output=", "lut-file="]

usage_string = """
Use this program to construct a simple or composite gas model for use with
the simulation codes Eilmer3 of L1d3.

Usage: 
  gasmodel.py [--help] [--model=<modelName>] [--species=<speciesList|none>] \\
              [--lut-file=<LUTFileName>] \\
              [--output=<luaFile|gas-model.lua>]

Input parameters:
  model   : name of the gas model, may have embedded spaces.
  species : list of species names (space delimited) in a single string.
  output  : name of the gas-model file to be written.
  lut-file: name of the preexisting LUT-gas model file, if relevant.

Examples:
  $ gasmodel.py --model='thermally perfect gas' --species='N2 N' \\
                --output='nitrogen.lua'

  $ gasmodel.py --model='ideal gas' --species='Ar He' \\
                --lut-file='cea-lut-custom.lua.gz' \\
                --output='LUT-plus-Ar-He.lua'

Notes:
  If you want a LUT-plus-composite gas model, set up the LUT table
  externally. Invoke this program, specifying the rest of the species
  for the composite gas model.  The LUT gas species is prepended to 
  the composite gas species list.  You will need both gas model files
  in place to use the resulting LUT-plus-composite gas model.
"""

def main():
    #
    try:
        userOptions = getopt(sys.argv[1:], [], longOptions)
    except GetoptError, e:
        print "One (or more) of your command-line options was no good."
        print "    ", e
        print usage_string
        sys.exit(1)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        print usage_string
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
