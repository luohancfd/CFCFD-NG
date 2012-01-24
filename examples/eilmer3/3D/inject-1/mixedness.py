#!/usr/bin/env python
# mixedness.py
#
# Custom postprocessor to compute a measure of mixedness plus a few other things
# (such as the average height, z_bar, of the fuel jet).
#
# PJ, 21-May-2010, 26-May-2010 (fixed enthalpy calc. and mass flux label)

import sys, os, gzip, math
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_grid import StructuredGrid
from e3_flow import StructuredGridFlow
from libprep3 import Vector, quad_centroid, quad_normal, quad_area, cross, dot

if len(sys.argv) != 3:
    print "Usage: mixedness.py jobName outputFileName"
    sys.exit()
jobName = sys.argv[1]
outputFileName = sys.argv[2]
tindx = '9999'

# We'll access the blocks using [i][j] indexing, with i progressing
# in the downstream direction.
nib = 3
njb = 2
bindx_list = [[0,1],[2,3],[4,5]]
flow = []
grd = []
for ib in range(nib):
    flow.append([])
    grd.append([])
    for jb in range(njb):
        bindx = '%04d' % bindx_list[ib][jb]
        print "bindx=", bindx
        fileName = 'flow/t'+tindx+'/'+ jobName + '.flow.b'+bindx+'.t'+tindx+'.gz'
        fp = gzip.open(fileName, "r")
        f = StructuredGridFlow()
        flow[-1].append(f)
        f.read(fp)
        fp.close()
        print "flow data: ni,nj,nk=", f.ni, f.nj, f.nk
        # The grid is always at tindx 0.
        fileName = 'grid/t0000/'+ jobName + '.grid.b'+bindx+'.t0000.gz'
        fp = gzip.open(fileName, "r")
        g = StructuredGrid()
        grd[-1].append(g)
        g.read(fp)
        fp.close()
        print "grid data: ni,nj,nk=", g.ni, g.nj, g.nk

# Work down the duct (i-direction) and compute fraction of area where gas is not
# purely one species or the other.
def process_one_slice(flow, grd, i):
    """
    Returns ...

    Use cell-average flow quantities and areas from the EAST boundary faces.
    
    This is a little dodgy, in general, but we gat away with it here
    because the cells have right angles.
    """
    TOL = 5.0e-2
    rho_f1_z_dA = 0.0; rho_z_dA = 0.0
    m_flux = 0.0; H_flux = 0.0
    A_f0 = 0.0; A_f1 = 0.0; A_mix = 0.0; A_tot = 0.0
    nj = flow.nj; nk = flow.nk
    for j in range(nj):
        for k in range(nk):
            vtx = grd.get_vertex_list_for_cell(i, j, k)
            # The EAST cell face has 1, 2, 6, 5 as corners (unit normal out).
            face_centroid = quad_centroid(vtx[1], vtx[2], vtx[6], vtx[5])
            face_normal = quad_normal(vtx[1], vtx[2], vtx[6], vtx[5])
            face_area = quad_area(vtx[1], vtx[2], vtx[6], vtx[5])
            # average conditions in cell
            z = flow.data['pos.z'][i,j,k]
            pressure = flow.data['p'][i,j,k] 
            rho = flow.data['rho'][i,j,k]
            e = flow.data['e[0]'][i,j,k]
            f0 = flow.data['massf[0]'][i,j,k]
            f1 = flow.data['massf[1]'][i,j,k]
            rho_z_dA += rho * face_area * z
            rho_f1_z_dA += rho * f1 * face_area * z
            A_tot += face_area
            if f0 > 1.0 - TOL:
                A_f0 += face_area
            elif f0 < TOL:
                A_f1 += face_area
            else:
                A_mix += face_area
            v_x = flow.data['vel.x'][i,j,k]
            v_y = flow.data['vel.y'][i,j,k]
            v_z = flow.data['vel.z'][i,j,k]
            vel_abs = Vector(v_x, v_y, v_z)
            df = face_area * pressure * face_normal
            dm_flux = rho * dot(vel_abs, face_normal) * face_area
            m_flux += dm_flux
            H_flux += dm_flux * (e + pressure/rho + 0.5*(v_x*v_x + v_y*v_y + v_z*v_z))
    return m_flux, H_flux, A_f0, A_f1, A_mix, A_tot, rho_z_dA, rho_f1_z_dA

fp = open(outputFileName, "w")
fp.write("# 1:x(m) 2:mass_flux(kg/s) 3:total_enthalpy_flux(W) ")
fp.write("4:A_f1(m**2) 5:A_f1(m**2) 6:A_mix(m**2) 7:A_tot(m**2) 8:A_mix/A_tot ")
fp.write("9:rho_z_dA 10:rho_f1_z_dA 11:z_bar(m)\n")

for ib in range(nib):
    for i in range(flow[ib][jb].ni):
        rho_f1_z_dA = 0.0; rho_z_dA = 0.0
        m_flux = 0.0; H_flux = 0.0
        A_f0 = 0.0; A_f1 = 0.0; A_mix = 0.0; A_tot = 0.0
        x = flow[ib][jb].data['pos.x'][i][0][0] # Assuming planar slices
        print x,
        sys.stdout.flush()
        for jb in range(njb):
            m, H, A0, A1, Amix, Atot, rhozdA, rhof1zdA = \
                process_one_slice(flow[ib][jb], grd[ib][jb], i)
            rho_z_dA += rhozdA; rho_f1_z_dA += rhof1zdA
            m_flux += m; H_flux += H
            A_f0 += A0; A_f1 += A1; A_mix += Amix; A_tot += Atot
        z_bar = rho_f1_z_dA/rho_z_dA
        fp.write("%g %g %g %g %g %g %g %g %g %g %g\n" % 
                 (x, m_flux, H_flux, A_f0, A_f1, A_mix, A_tot, A_mix/A_tot,
                  rho_z_dA, rho_f1_z_dA, z_bar))

fp.close()
print "Done."


