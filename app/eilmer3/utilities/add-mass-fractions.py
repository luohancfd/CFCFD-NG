#!/usr/bin/env python

import sys
from numpy import zeros
from e3_flow import read_all_blocks, write_Tecplot_file
from cfpylib.gasdyn.cea2_gas import Gas

nblocks = 32
speciesList = ['N2', 'O2', 'NO', 'N', 'O']

def print_usage():
    print "Usage: ./add_mass_fractions.py jobName tindx"
    print ""
    print "where:"
    print "  jobName -- base name for job"
    print "  tindx   -- time index of flow solution"
    print ""

def read_time_indices(rootName):
    """
    The job.times should contain the map from solution index to simulation time.

    The file format is one entry per line, possibly with comment lines 
    starting with the # character.

    Returns the final tindx and a dictionary of (tindx,t)-pairs
    """
    fileName = rootName + ".times"
    try:
        fp = open(fileName, "r")
        # For unkown reason, sometimes a line has several
        # leading null-characters. So remove them to avoid
        # an error while reading the file. (Stefan Hess)
        buf = fp.readline().strip().strip("\0")
        times_dict = {}
        while len(buf) > 0:
            if buf[0] == '#': 
                buf = fp.readline().strip().strip("\0")
                continue
            tokens = buf.split()
            tindx = int(tokens[0]); t = float(tokens[1])
            times_dict[tindx] = t
            buf = fp.readline().strip().strip("\0")
        fp.close()
    except Exception, e:
        # For some reason, we could not read the job.times file
        # so let's assume that there should be a t=0 data set.
        # This will allow us to handle cases where we have run
        # the preparation program but not the main simulation.
        print "Problem reading .times file:", e
        tindx = 0; times_dict = { 0:0.0 }
    return tindx, times_dict


def add_mass_fractions_to_block(blk):
    # Create cea2_gas object
    air = Gas({'Air':1.0}, outputUnits='massf')
    for isp, s in enumerate(speciesList):
        key = "massf[%i]-%s" % (isp, s)
        blk.data[key] = zeros((blk.ni,blk.nj,blk.nk),'d')
        blk.vars.append(key)

    for k in range(blk.nk):
        for j in range(blk.nj):
            for i in range(blk.ni):
                p = blk.data['p'][i,j,k]
                T = blk.data['T[0]'][i,j,k]
                air.set_pT(p, T)
                massf = air.get_fractions(speciesList)
                for isp, mf in enumerate(massf):
                    key = "massf[%i]-%s" % (isp, speciesList[isp])
                    blk.data[key][i,j,k] = mf
    return

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print_usage()
        sys.exit(1)

    print "add-mass-fractions : Adding mass fractions to equilibrium air simulation."
    
    job = sys.argv[1]
    tindx = int(sys.argv[2])
    final_tindx, times_dict = read_time_indices(job)

    print "Reading flow data."
    grid, flow, dimensions = read_all_blocks(job, nblocks, tindx, True, False)
    print "Done."

    print "Adding mass fractions using CEA."
    for blk in flow:
        add_mass_fractions_to_block(blk)
    print "Done."
    
    print "Writing out Tecplot file."
    write_Tecplot_file(job, tindx, nblocks, grid, flow, times_dict[tindx])
    print "Done."

    print "add-mass-fractions : Done."

        
    
        
