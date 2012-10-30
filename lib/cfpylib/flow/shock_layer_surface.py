"""
shock_layer_surface.py: A class for generating a 2D surface between two
                        curves ensuring that the cells adjacent to the 
                        east curve are normal to the surface.

Author: DFP

Version: 30-Oct-2012 : initial repository version, ported from my local version.
"""


try:
    from libprep3 import *
except:
    print "Could not import libprep3."
    sys.exit()
from math import *
from cfpylib.nm.zero_solvers import bisection

class ShockLayerSurface:
    def __init__( self, east, west ):
	self.east = east
	self.west = west

    def eval( self, r, s ):
	# eval the location on the east face
	ep = self.east.eval(s)
	# eval normal angle on the east face
	dpdt = self.east.dpdt(s)
	if dpdt.x == 0.0:
            theta = 0.0
        else:
  	    theta = - atan( dpdt.y / dpdt.x ) + pi / 2.	
	#  find the west point
	if s==0.0: t = 0
        elif s==1.0: t = 1
	elif theta != 0.0:
            def f( t ):
	        wp = self.west.eval(t)
	        L = ( wp.y - ep.y ) / sin( theta ) 
	        lp = Vector3( ep.x - L * cos(theta), ep.y + L * sin(theta) ) 
	        return wp.x - lp.x
	    t = bisection( f, 0.0, 1.0 )
	else:
            def f( t ):
	        wp = self.west.eval(t)
	        L = ( ep.x - wp.x ) / cos( theta ) 
	        lp = Vector3( ep.x - L * cos(theta), ep.y + L * sin(theta) ) 
	        return wp.y - lp.y
	    t = bisection( f, 0.0, 1.0 )
	wp = self.west.eval(t)
        # make the normal line
	line = Line(ep,wp)
	# eval the point
	return line.eval(1-r)

    def str():
	return "ShockLayerSurface"

