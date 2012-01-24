#! /usr/bin/env python
# \file estimate_ramp_force.py
#
# Example postprocessing script to look at the data along the ramp
# and compute some potentially useful information.

import sys, os, string
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in working directory
from e3_grid import StructuredGrid
from e3_flow import StructuredGridFlow
from libprep3 import *
from gzip import GzipFile

print "\n\nEstimate force on the ramp surface."

fileName = 'grid/t0000/simple_ramp.grid.b0001.t0000.gz'
print "Read grid file:", fileName
fin = GzipFile(fileName, "rb")
grd = StructuredGrid()
grd.read(f=fin)
fin.close()
print "Read grid: ni=", grd.ni, "nj=", grd.nj, "nk=", grd.nk

fileName = 'flow/t0005/simple_ramp.flow.b0001.t0005.gz'
print "Read solution file:", fileName
fin = GzipFile(fileName, "rb")
soln = StructuredGridFlow()
soln.read(fin)
fin.close()
ni = soln.ni; nj = soln.nj; nk = soln.nk
print "Read solution: ni=", ni, "nj=", nj, "nk=", nk

# Integrate the pressure force over the BOTTOM surface of the block.
force = Vector(0.0, 0.0, 0.0)
k = 0  
for i in range(ni):
    for j in range(nj):
        p0,p1,p2,p3,p4,p5,p6,p7 = grd.get_vertex_list_for_cell(i,j,k)
        # The bottom cell face has p0, p1, p2, p3 as corners.
        surface_centroid = quad_centroid(p0, p1, p2, p3)
        surface_normal = quad_normal(p0, p1, p2, p3)
        surface_area = quad_area(p0, p1, p2, p3)
        pressure = soln.data["p"][i][j][k] # average pressure in cell
        df = surface_area * pressure * surface_normal
        force -= df # negative because the unit normal of this cell face is into the volume
print "force=", force, "Newtons"

# Find the distance from the cell centre to the centroid of the cell face
# for a strip of cells along the ramp.  Although it is not of much use here,
# this information could be used to estimate the boundary-layer growth
# along the plate.  Katsu did this for his scramjet calculations.
fileName = "distances.txt"
fout = open(fileName, "w")
k = 0; j = 0
for i in range(ni):
    p0,p1,p2,p3,p4,p5,p6,p7 = grd.get_vertex_list_for_cell(i,j,k)
    # The bottom cell face has p0, p1, p2, p3 as corners.
    surface_centroid = quad_centroid(p0, p1, p2, p3)
    surface_normal = quad_normal(p0, p1, p2, p3)
    surface_area = quad_area(p0, p1, p2, p3)
    # We pull the cell-centre information out of the solution data.
    cell_centre = Vector(soln.data["pos.x"][i][j][k],
                         soln.data["pos.y"][i][j][k],
                         soln.data["pos.z"][i][j][k])
    distance1 = vabs(cell_centre - surface_centroid)
    # We compute the cell centroid from the grid vertices.
    cell_centre_2 = hexahedron_centroid(p0, p1, p2, p3, p4, p5, p6, p7)
    distance2 = vabs(cell_centre_2 - surface_centroid)
    fout.write("%d %e %e\n" % (i, distance1, distance2))
fout.close()

print "done."
