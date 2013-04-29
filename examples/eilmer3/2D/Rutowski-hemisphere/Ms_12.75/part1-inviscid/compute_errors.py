#!/usr/bin/env python

from cfpylib.util.YvX import YvX
from numpy import array, sqrt

def get_rms( y_list ):
    y_array = array(y_list)
    return sqrt(sum(y_array**2)/len(y_array))
    
print "\nPerforming comparison with solution of Marco Panesi (VKI PhD thesis 2009)...\n"

# ref data
ref_TvX = YvX( "ref_profile.data", 0, 24 )

# simulation data
sim_TvX = YvX( "profile.data", 0, 24 )

# error in the vibration-electron-electronic temperature
error_in_T = []
for i,x in enumerate(sim_TvX.x_array[:-2]):
    T_ref = ref_TvX.y_from_x(x)
    T_sim = sim_TvX.y_array[i]
    error_in_T.append( ( T_sim - T_ref ) / T_ref )
    
print "RMS of the T_ve error is %0.2f percent" % ( get_rms( error_in_T )*100 )
