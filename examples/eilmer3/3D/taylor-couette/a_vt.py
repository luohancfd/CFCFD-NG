#! /usr/bin/env python
# \file a_vt.py
#
# caculate the averaged velocity and temperature along the radial gap

import sys, os, string
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in working directory
from e3_grid import StructuredGrid
from e3_flow import StructuredGridFlow
from libprep3 import *
from gzip import GzipFile
from math import sin, cos, tan, atan, pi, sqrt

print "\n\ncaculate the average temperature and velocity."

fileName = 'grid/t0000/tc_flow_nitrogen.grid.b0000.t0000.gz'
print "Read grid file:", fileName
fin = GzipFile(fileName, "rb")
grd = StructuredGrid()
grd.read(f=fin)
fin.close()
print "Read grid: ni=", grd.ni, "nj=", grd.nj, "nk=", grd.nk

fileName = 'flow/t0036/tc_flow_nitrogen.flow.b0000.t0036.gz'
print "Read solution file:", fileName
fin = GzipFile(fileName, "rb")
soln = StructuredGridFlow()
soln.read(fin)
fin.close()
ni = soln.ni; nj = soln.nj; nk = soln.nk
print "Read solution: ni=", ni, "nj=", nj, "nk=", nk

# Caculate the averaged velocity and temperature along the radial gap
# South surface of block 0
fileName = "average.txt"
fout = open(fileName, "w")
j = 0  
v_tan = 0.0
vel_tan = 0.0
T_tan = 0.0
Tem_tan = 0.0
for i in range(ni):
    for k in range(nk):
        pos_x = soln.data["pos.x"][i][j][k]
        pos_y = soln.data["pos.y"][i][j][k] 
        r_g = sqrt(pos_x*pos_x + pos_y*pos_y)
        r_1 = ( r_g-0.2125 ) / 0.0031                # radial position
        vel_x = soln.data["vel.x"][i][j][k]
        vel_y = soln.data["vel.y"][i][j][k]
        T_tan = soln.data["T[0]"][i][j][k]           # temp temperature
        theta = atan(pos_y/pos_x)
        v_tan = vel_y*cos(theta) - vel_x*sin(theta)  # temp tangential velocity
        v_tan = v_tan / 614.0                        # change into nondimensionalized form
        vel_tan += v_tan                             # the sum of tangential velocity
        Tem_tan += T_tan                             # the sum of tempearture
    fout.write("%e %e %e\n" % (r_1, vel_tan/nk, Tem_tan/nk))
    v_tan = 0.0
    vel_tan = 0.0
    T_tan = 0.0
    Tem_tan = 0.0
fout.close()

print "done."
