#!/usr/bin/env python

from cfpylib.util.YvX import YvX
from numpy import array, sqrt

def get_rms( y_list ):
    y_array = array(y_list)
    return sqrt(sum(y_array**2)/len(y_array))
    
print "\nPerforming comparison with solution of Marco Panesi (VKI PhD thesis 2009)...\n"

# UTIAS experimental data (ionization fraction) from Glass and Liu
exp_AvX = YvX( "glass_liu_ne_M_16.5.txt" )

# my data (electron mole fraction = ionization fraction)
p3_AvX = YvX( "output.data", 0, 8 )
p3_AvX.x_array[:] *= 100.0
p3_AvX.recompute_spline()

# error in the ionization fraction
error_in_alpha = []
for i,x in enumerate(p3_AvX.x_array):
    alpha_exp = exp_AvX.y_from_x(x)
    alpha_p3 = p3_AvX.y_array[i]
    error_in_alpha.append( ( alpha_p3 - alpha_exp ) / alpaha_exp )
    
print "RMS of the alpha error is %0.2f percent" % ( get_rms( error_in_alpha )*100 )
