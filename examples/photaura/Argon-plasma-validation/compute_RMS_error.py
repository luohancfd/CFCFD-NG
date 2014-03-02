#!/usr/bin/env python
# compute_RMS_error.py

import sys, os
from cfpylib.util.YvX import YvX
import math

def printUsage():
    print "compute_RMS_error.py"
    print "Compute the RMS error between two datasets."
    print "Usage:"
    print "compute_RMS_error.py FILE1 XCOL1 YCOL1 FILE2 XCOL2 YCOL2"
    sys.exit(1)
       
def main():
    if len(sys.argv)!=7:
	printUsage()
    
    sol = YvX( sys.argv[1], int(sys.argv[2])-1, int(sys.argv[3])-1 ) 
    ref = YvX( sys.argv[4], int(sys.argv[5])-1, int(sys.argv[6])-1 )
    
    y_RMS = 0.0
    count = 0
    for i in range(len(sol.x_array)):
        x = sol.x_array[i]
        y = sol.y_array[i]
	y_ref = ref.y_from_x(y)
        error = ( y - y_ref ) / y_ref
	print error
	y_RMS += error**2
        count  += 1
    
    y_RMS = math.sqrt( y_RMS / count )
            
    print "RMS error = %0.2f percent" % ( y_RMS * 100 )

    print "Done.\n"

if __name__ == '__main__':
    main()
