#! /usr/bin/env python
##
## \file wall_definition.py
## \ingroup geom
##
## \brief Define a nozzle wall in terms of nodes of a Bezier polynomial.
##
## \author PA Jacobs
##
## \version 21-Mar-05 first attempt, X2 nozzle

from Numeric import array, arange
from roberts import distribute_points_1
from bezier import bezier3D_eval

def bezier_nodes_for_wall(parameter_list):
    """Generates the coordinates for the nodes defining the wall."""
    y_inlet = 0.0425  # inner radius of X2 acceleration tube
    n_nodes = 1 + len(parameter_list) + 1
    xmax = 1.2  # metres
    beta = 1.1  # clustering parameter (see Roberts stretching function)
    x = distribute_points_1(0.0, xmax, n_nodes-1, 1, 0, beta)
    y = array([y_inlet,]* n_nodes)
    for i in range(1,n_nodes):
        iparam = i - 1
        if iparam < len(parameter_list):
            y[i] = parameter_list[iparam]
        else:
            # the last y-node is at the same height as the
            # second-last y-node to get a zero-slope exit.
            y[i] = parameter_list[-1]
    return x, y

#------------------------------------------------------------------------

if __name__ == '__main__':
    print "# Begin demo of wall_definition for X2 nozzle..."
    p_list= [0.040496, 0.065839, 0.051086, 0.089875, 0.103871]
    x, y = bezier_nodes_for_wall(p_list)
    print "#---------------------"
    print "# Bezier nodes..."
    for i in range(len(x)):
        print x[i], y[i]
    print "#---------------------"
    print "# Sample points..."
    for t in arange(0.0,1.001,0.02):
        xp, yp, zp = bezier3D_eval(t, Bx=x, By=y)
        print xp, yp
    print "#---------------------"
    print "Done."
