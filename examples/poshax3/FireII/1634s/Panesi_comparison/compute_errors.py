#!/usr/bin/env python

from cfpylib.util.YvX import YvX
from numpy import array, sqrt

def get_rms( y_list ):
    y_array = array(y_list)
    return sqrt(sum(y_array**2)/len(y_array))
    
print "\nPerforming comparison with solution of Marco Panesi (VKI PhD thesis 2009)...\n"

# panesi data
panesi_TvX = YvX( "marco_1634s_Tv.txt" )

# my data
p3_TvX = YvX( "output.data", 0, 2 )

# error in the vibration-electron-electronic temperature
error_in_T = []
for i,x in enumerate(p3_TvX.x_array[:-2]):
    T_panesi = panesi_TvX.y_from_x(x)
    T_p3 = p3_TvX.y_array[i]
    error_in_T.append( ( T_p3 - T_panesi ) / T_panesi )
    
print "RMS of the T_ve error is %0.2f percent" % ( get_rms( error_in_T )*100 )
