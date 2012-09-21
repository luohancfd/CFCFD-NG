#!/usr/bin/env python
# nenzfr_compute_viscous_data.py
# 
# Luke Doherty
# 16-Aug-2011
#
# History:
#   This script is based directly on Wilson Chan's estimate-skin-friction.py 
#  script which was in turn based on Peter Jacob's compute_y_plus.py script. 
#   Wilson's script was re-written to make it suitable for use with 2D(nozzle)
#  geometries and then further modified for inclusion/integration with 
#  nenzfr.py.
#                                                         Luke D. 16-Aug-2011

VERSION_STRING = "Luke D, 27-May-2012"

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from e3_flow import *
from libprep3 import *
import optparse

#----------------------------------------------------------------------------
def get_cell_corners(grd, i, j, k):
    p0 = Vector(grd.x[i][j][k], grd.y[i][j][k], grd.z[i][j][k])
    p1 = Vector(grd.x[i+1][j][k], grd.y[i+1][j][k], grd.z[i+1][j][k])
    p2 = Vector(grd.x[i+1][j+1][k], grd.y[i+1][j+1][k], grd.z[i+1][j+1][k])
    p3 = Vector(grd.x[i][j+1][k], grd.y[i][j+1][k], grd.z[i][j+1][k])
    return p0, p1, p2, p3

#----------------------------------------------------------------------------

def compute_viscous(jobName, nblock, nbj, tindx):
    # Read in all grid and flow files
    zipFiles = 1
    grid, flow, bgk, dimensions = read_all_blocks(jobName, nblock, tindx, zipFiles)
    
    for jb in range(nbj-1, nblock, nbj):
        if jb==nbj-1:
            outfile = open(jobName+"-viscous.data", "w")
            outfile.write("#x(m)\ty(m)\tz(m)\tdu(m/s)\tdy(m)\trho(kg/m^3)\tmu(N.s/m^2)\ttau_w(N/m^2)\ty_plus \n")
        else:
            outfile = open(jobName+"-viscous.data", "a")

        # For NORTH wall
        j = flow[jb].nj - 1
        for i in range(flow[jb].ni):
            for k in range(flow[jb].nk):
                x = flow[jb].data['pos.x'][i,j,k]
                y = flow[jb].data['pos.y'][i,j,k]
                z = flow[jb].data['pos.z'][i,j,k]
                
                vel_vector = Vector(flow[jb].data['vel.x'][i,j,k], flow[jb].data['vel.y'][i,j,k],\
                                    flow[jb].data['vel.z'][i,j,k])
                du = vabs(vel_vector)
                
                p0, p1, p2, p3 = get_cell_corners(grid[jb], i, j, k)            
                
                # The top face has p3, p2 as corners.
                edge_center = (p3+p2)/2
                edge_vector = p3-p2
                unit_normal = Vector(-edge_vector.y, edge_vector.x, 0)/vabs(edge_vector)
                    
                dy = abs(dot((Vector(x, y, z) - edge_center), unit_normal))
                rho = flow[jb].data['rho'][i,j,k]
                mu = flow[jb].data['mu'][i,j,k]
                tau_w = mu * du/dy   # Wall shear stress
                
                u_tau = (tau_w / rho)**0.5   # Friction velocity
                yplus = u_tau * dy * rho / mu
	        outfile.write("%1.6e\t%1.6e\t%f\t%1.6f\t%1.6e\t%1.6e\t%1.6e\t%0.3f\t%f \n" % (x, y, z, du, dy, rho, mu, tau_w, yplus))
        outfile.close()


def main():
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('--job', dest='jobName', default=None,
                  help="base name for Eilmer3 files [default: %default]")
    op.add_option('--tindx', dest='tindx', type='int', default=9999, 
                  help="time index for solution data to be processed [default: %default]")
    op.add_option('--nblock', dest='nblock', type='int', default=None,
                  help="total number of blocks to consider")
    op.add_option('--nbj', dest='nbj', type='int', default=1,
                  help="number of blocks in the radial/transverse direction "
                       "[default: %default]")
    opt, args = op.parse_args()
    
    # If no value for nblock is given, we infer the number of blocks from the number of 
    # solution files found in the './flow' directory.
    if opt.nblock is None:
        flowdirectory = os.getcwd() + "/flow/t" + str(opt.tindx) + "/"
        opt.nblock = len( os.listdir(flowdirectory) )

    # Check that a jobName has been given.
    bad_input = False
    if opt.jobName is None:
        print "Need to supply a jobName"
        bad_input = True
    if bad_input:
        return -2    
    
    # Compute y_plus data.
    compute_viscous(opt.jobName, opt.nblock, opt.nbj, opt.tindx)
    
    print "Computed viscous data may be found in "+opt.jobName+"-viscous.data"
    return 0


#----------------------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr-compute-viscous-data"
        print "   Version:", VERSION_STRING
        print "   For more details use option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
