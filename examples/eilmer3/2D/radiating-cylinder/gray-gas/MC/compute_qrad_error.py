#!/usr/bin/env python
# compute_qrad.py

import sys, os
from cfpylib.util.YvX import YvX

def printUsage():
    print "compute_qrad.py"
    print "Compute the average radiative heat flux from a file between two spatial limits."
    print "Usage:"
    print "compute_qrad_error.py FILE1 XCOL1 YCOL1 S1 S2 QRAD"
    sys.exit(1)
       
def main():
    if len(sys.argv)!=7:
	printUsage()
    
    sol = YvX( sys.argv[1], int(sys.argv[2])-1, int(sys.argv[3])-1 ) 
    s1 = float(sys.argv[4])
    s2 = float(sys.argv[5])
    
    
    s_list = []
    q_list = []
    for i in range(len(sol.x_array)):
        s = sol.x_array[i]
        q = sol.y_array[i]
        if s >= s1 and s <= s2:
            s_list.append(s)
            q_list.append(q)
            
    q_av = sum( q_list ) / len( q_list ) 
    
    q_rad = float( sys.argv[6] ) 
    error = ( q_av - q_rad ) / q_av
        
    print "qrad error = %0.2f percent" % ( error * 100 )

    print "Done.\n"

if __name__ == '__main__':
    main()
