# file: compute_shock_standoff.py
#
# AP, 28-Feb-2013

import numpy as np
from e3_flow import *
from libprep3 import *
                    
def get_boundary_position(grids, boundary=3):
    X = np.array([])
    Y = np.array([])
    for grid in grids:
        if boundary == 0: # North
            x = grid.x[:,-1,0]
            y = grid.y[:,-1,0]
        elif boundary == 1: # East
            x = grid.x[-1,:,0]
            y = grid.y[-1,:,0]
        elif boundary == 2: # South
            x = grid.x[:,0,0]
            y = grid.y[:,0,0]
        elif boundary == 3: # West
            x = grid.x[0,:,0]
            y = grid.y[0,:,0]
        else:
            print "Boundary %s not defined." % (boundary)
            x = []
            y = []
            
        X = np.append(X,x)
        Y = np.append(Y,y)
    return X, Y
    
rootName = 'sphere3'
grids, g, h, j = read_all_blocks(rootName, 4, 9999, zipFiles=1, movingGrid=1)
x3,y3 = get_boundary_position(grids) # West boundary, i.e. the shock
x1,y1 = get_boundary_position(grids, boundary=1) # East boundary, i.e. the body

print 'shock_standoff_distance=', np.abs(x3[0]-x1[0])


