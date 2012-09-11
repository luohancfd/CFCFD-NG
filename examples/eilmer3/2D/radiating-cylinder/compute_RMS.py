#!/usr/bin/env python
# compute_relative_errors.py

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in working directory
from cfpylib.util.YvX import YvX
from math import sqrt

make_plots = False

def printUsage():
    print "compute_RMS.py"
    print "Compute the relative errors between columns of values from two different files."
    print "NOTE: first file is the data of interest, second file is the reference solution."
    print "Usage:"
    print "compute_relative_errors.py FILE1 XCOL1 YCOL1 FILE2 XCOL2 YCOL2"
    sys.exit(1)
   
def compute_RMS_error( sol, ref ):
    tmp = 0.0
    n = len(sol.x_array)
    for i in range(n):
        x = sol.x_array[i]
        y_ref = ref.y_from_x(x)
        dy = sol.y_array[i] - y_ref
        tmp += (dy/y_ref)**2
	
	return sqrt( tmp / n )
       
def main():
    if len(sys.argv)!=7:
	printUsage()
    
    sol = YvX( sys.argv[1], int(sys.argv[2])-1, int(sys.argv[3])-1 ) 
    ref = YvX( sys.argv[4], int(sys.argv[5])-1, int(sys.argv[6])-1 )
    
    rms_error = compute_RMS_error( sol, ref )
    
    print "RMS error = %0.2f percent" % ( rms_error * 100 )
    
    if make_plots:
        sol.plot_data( title="radiative divergence across center of slab", xlabel="distance", ylabel="radiative divergence", label="e3: MC", new_plot=True, show_plot=False, include_integral=False, rep='-', logscale_y=False, xrange=None, yrange=None ) 
        ref.plot_data( title="radiative divergence across center of slab", xlabel="distance", ylabel="radiative divergence", label="exact", new_plot=False, show_plot=True, include_integral=False, rep='-', logscale_y=False, xrange=None, yrange=None )        
    
    print "Done.\n"

if __name__ == '__main__':
    main()
