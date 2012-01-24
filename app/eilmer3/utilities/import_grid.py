#! /usr/bin/env python
## \file import_grid.py
## \ingroup elmer
## \brief Import grid specification from a given format into Eilmer3 format.
##
## \author PJ
##
## \version 02-Dec-2004 all in one file
## \version 05-Apr-2005 most of the smarts for a single block grid
##                      have been moved to the separate module grid.py
## \version 07-Feb-2010 port to Eilmer3
##
    
from getopt import getopt
import sys, os, string
import sys
import os
sys.path.append("/sw/lib/python2.3/site-packages/Numeric") # for Tim's MacOSX
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory

from e3_grid import StructuredGrid

shortOptions = ""
longOptions = ["input=", "output=", "singleblock",
               "plot3dplanes", "plot3dwhole", "noblanking"]

def printUsage():
    print ""
    print "Usage: import_grid.py" + \
          " --input=<inputFileName>" + \
          " --output=<baseName>" + \
          " [--plot3dplanes|--plot3dwhole]" + \
          " [--singleblock]" + \
          " [--noblanking]"
    print ""
    return
    
if __name__ == '__main__':
    print "Begin import_grid..."
    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) < 2:
        printUsage()
        sys.exit()
        
    inputFileName = uoDict.get("--input", "")
    print "Read grid from", inputFileName
    baseFileName = uoDict.get("--output", "test")
    print "Construct block file names from", baseFileName
    
    fin = open(inputFileName, "r")

    # For multiple-block plot3D files,
    # the first line contains the number of blocks. 
    if uoDict.has_key("--singleblock"):
        nblock = 1
    else:
        lineContent = fin.readline()
        nblock = int(lineContent)
    print "Expected number of blocks: ", nblock

    # There are two flavours: data in planes or whole-block
    dataInPlanes = 0;  # default
    if uoDict.has_key("--plot3dplanes"):
        dataInPlanes = 1
    if uoDict.has_key("--plot3dwhole"):
        dataInPlanes = 0

    # Some grid-generation software includes blanking data, some doesn't.
    if uoDict.has_key("--noblanking"):
        with_blanking = 0
    else:
        with_blanking = 1
        
    print "Collect numbers of indices"
    indices = []
    while len(indices) < nblock * 3:
        lineContent = fin.readline()
        words = lineContent.split()
        for word in words:
            indices.append(int(word))

    print "Start reading the blocks."
    blocks = []
    for ib in range(nblock):
        print "Block:", ib
        ni = indices[ib * 3]
        nj = indices[ib * 3 + 1]
        nk = indices[ib * 3 + 2]
        b = StructuredGrid((ni, nj, nk))
        if dataInPlanes:
            b.read_from_plot3d_in_planes(fin, with_blanking)
        else:
            b.read_from_plot3d_whole_grid(fin, with_blanking)
        vtkFileName = baseFileName + "." + str(ib) + ".g.vtk"
        b.label = baseFileName + ": block=" + str(ib)
        fout = open(vtkFileName, "w")
        b.write_block_in_VTK_format(fout)
        fout.close()
    
    print "Done."
