#!/usr/bin/env python
# get_residuals.py

import sys, os

def printUsage():
    print "get_residuals.py"
    print "Pull out residuals from log files and make list against time (and steps) in a new file"
    print "Usage:"
    print "get_residuals.py BLOCK OUTFILE"
    print "BLOCK: either block index or 'global'"
    sys.exit(1)
    
class ResidualVsTime:
    def __init__( self, s, t, rr, re ):
        self.s = s
        self.t = t
        self.rr = rr
        self.re = re
    def str(self):
        ostr = ""
        for i in range(len(self.s)):
            ostr += "%d \t %e \t %e \t %e\n" % ( self.s[i], self.t[i], self.rr[i], self.re[i] )
        return ostr

def main():
    if len(sys.argv)<3: printUsage()
    
    if sys.argv[1]=="global":
        ib = 0
    else:
        ib = int(sys.argv[1])
    block_tag = sys.argv[1]

    logfile_name = "e3mpi.%04d.log" % ib
    print logfile_name
    if ( os.path.isfile(logfile_name) ):
        print "Found log file for %d block" % ib   
        inFile = open( logfile_name, "r" )
        steps = []; times = []; density_residuals = []; energy_residuals = []
        while 1:
            line = inFile.readline()
            if len(line)==0: break
            tks = line.split()
            if block_tag!="global":
                if tks[0]=="RESIDUAL" and tks[1]=="mass" and tks[2]=="block":
                    density_residuals.append( float(tks[5]) )
                if tks[0]=="RESIDUAL" and tks[1]=="energy" and tks[2]=="block":
                    energy_residuals.append( float(tks[5]) )
                if tks[0]=="RESIDUAL" and tks[1]=="mass" and tks[2]=="global":
                    steps.append(int(tks[6]))
                    times.append(float(tks[8]))
            else:
                if tks[0]=="RESIDUAL" and tks[2]=="global":
                    if tks[1]!="mass":
                        print line
                        print "expected this to be the density residual line!"
                        sys.exit()
                    steps.append( int(tks[6]) )
                    times.append( float(tks[8]) )
                    density_residuals.append(float(tks[4]))
                    line = inFile.readline()
                    tks = line.split()
                    if tks[1]!="energy":
                        print line
                        print "expected this to be the energy residual line!"
                        sys.exit()
                    energy_residuals.append(float(tks[4]))
        residuals = ResidualVsTime( s=steps, t=times, rr=density_residuals, re=energy_residuals )
        inFile.close()
    else:
        print "Could not find requested logfile: ", logfile_name
        sys.exit()
	
    outfile_name = sys.argv[2]
    outfile = open(outfile_name, 'w')
    outfile.write("# This file was created by get_residuals.py\n")
    outfile.write("# Residuals for block: %s\n" % block_tag )
    outfile.write("# Column 1: Step, s (ND)\n")
    outfile.write("# Column 2: Time, t (seconds)\n")
    outfile.write("# Column 3: Density residual (ND)\n")
    outfile.write("# Column 4: Energy residual (ND)\n")
    outfile.write(residuals.str())

    outfile.close()
    
    print "Done.\n"

if __name__ == '__main__':
    main()
