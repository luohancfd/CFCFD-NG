#!/usr/bin/env python
# mix_and_vort.py
#
# Custom postprocessor to compute a measure of mixedness plus a few other things
# (such as jet penetration, total pressure loss and streamwise circulation).
#
# PJ, 21-May-2010 & VW, 26-July-2010

import sys, os, gzip, math
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_grid import StructuredGrid
from e3_flow import StructuredGridFlow
from libprep3 import Vector, quad_centroid, quad_normal, quad_area, cross, dot

if len(sys.argv) != 3:
    print "Usage: mix_and_vort.py jobName outputFileName"
    sys.exit()
jobName = sys.argv[1]
outputFileName = sys.argv[2]
tindx = '9999'
# tindx = '0014'

# We'll access the blocks using [i][j] indexing, with i progressing
# in the downstream direction.
nib = 6
njb = 4
bindx_list = [[0,1,2,3],[4,6,7,8],[5,9,10,11],[12,15,16,17],[13,18,19,20],[14,21,22,23]]
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
    """
    m_flux = 0.0; mf_flux=0.0; mfR_flux=0.0; H_flux = 0.0
    Circ = 0.0; prec = 0.0; z_pen = 0.0;
    nj = flow.nj; nk = flow.nk
    f0stoic = 0.02876
    f0lim = 0.1*f0stoic	# limiting mass fraction defining edge of jet
    for j in range(nj):
        for k in range(nk):
            vtx = grd.get_vertex_list_for_cell(i, j, k)
            # The EAST cell face has 1, 2, 6, 5 as corners (unit normal out).
            face_centroid = quad_centroid(vtx[1], vtx[2], vtx[6], vtx[5])
            face_normal = quad_normal(vtx[1], vtx[2], vtx[6], vtx[5])
            face_area = quad_area(vtx[1], vtx[2], vtx[6], vtx[5])
            # average conditions in cell
            y = flow.data['pos.y'][i,j,k]
            z = flow.data['pos.z'][i,j,k]
            pressure = flow.data['p'][i,j,k] 
            rho = flow.data['rho'][i,j,k]
	    e = flow.data['e[0]'][i,j,k]
            a = flow.data['a'][i,j,k]
	    v_x = flow.data['vel.x'][i,j,k]
            v_y = flow.data['vel.y'][i,j,k]
            v_z = flow.data['vel.z'][i,j,k]
	    if j==0:
	        if k==0:
                    yp = flow.data['pos.y'][i,j+1,k]
                    zp = flow.data['pos.z'][i,j,k+1]
                    v_yp = flow.data['vel.y'][i,j,k+1]
                    v_zp = flow.data['vel.z'][i,j+1,k]
	            vort_x = (v_zp-v_z)/(yp-y) - (v_yp-v_y)/(zp-z)
		elif k==nk-1:
                    yp = flow.data['pos.y'][i,j+1,k]
                    zm = flow.data['pos.z'][i,j,k-1]
                    v_zp = flow.data['vel.z'][i,j+1,k]
                    v_ym = flow.data['vel.y'][i,j,k-1]
	            vort_x = (v_zp-v_z)/(yp-y) - (v_y-v_ym)/(z-zm)
		else:
                    yp = flow.data['pos.y'][i,j+1,k]
                    zp = flow.data['pos.z'][i,j,k+1]
                    zm = flow.data['pos.z'][i,j,k-1]
                    v_zp = flow.data['vel.z'][i,j+1,k]
                    v_yp = flow.data['vel.y'][i,j,k+1]
                    v_ym = flow.data['vel.y'][i,j,k-1]
	            vort_x = 0.5 * ( 2.0*(v_zp-v_z)/(yp-y) - (v_yp-v_y)/(zp-z) - (v_y-v_ym)/(z-zm) )
	    elif j==nj-1:   
	        if k==0:
                    ym = flow.data['pos.y'][i,j-1,k]
                    zp = flow.data['pos.z'][i,j,k+1]
                    v_yp = flow.data['vel.y'][i,j,k+1]
                    v_zm = flow.data['vel.z'][i,j-1,k]
	            vort_x = (v_zm-v_z)/(ym-y) - (v_yp-v_y)/(zp-z)
		elif k==nk-1:
                    ym = flow.data['pos.y'][i,j-1,k]
                    zm = flow.data['pos.z'][i,j,k-1]
                    v_zm = flow.data['vel.z'][i,j-1,k]
                    v_ym = flow.data['vel.y'][i,j,k-1]
	            vort_x = (v_zm-v_z)/(ym-y) - (v_y-v_ym)/(z-zm)
		else:
                    ym = flow.data['pos.y'][i,j-1,k]
                    zp = flow.data['pos.z'][i,j,k+1]
                    zm = flow.data['pos.z'][i,j,k-1]
                    v_zm = flow.data['vel.z'][i,j-1,k]
                    v_yp = flow.data['vel.y'][i,j,k+1]
                    v_ym = flow.data['vel.y'][i,j,k-1]
	            vort_x = 0.5 * ( 2.0*(v_zm-v_z)/(ym-y) - (v_yp-v_y)/(zp-z) - (v_y-v_ym)/(z-zm) )
            elif k==0:
                yp = flow.data['pos.y'][i,j+1,k]
                zp = flow.data['pos.z'][i,j,k+1]
                ym = flow.data['pos.y'][i,j-1,k]
                v_yp = flow.data['vel.y'][i,j,k+1]
                v_zp = flow.data['vel.z'][i,j+1,k]
                v_zm = flow.data['vel.z'][i,j-1,k]
	        vort_x = 0.5 * ( (v_zp-v_z)/(yp-y) + (v_z-v_zm)/(y-ym) - 2.0*(v_yp-v_y)/(zp-z) )
	    elif k==nk-1:
                yp = flow.data['pos.y'][i,j+1,k]
                ym = flow.data['pos.y'][i,j-1,k]
                zm = flow.data['pos.z'][i,j,k-1]
                v_zp = flow.data['vel.z'][i,j+1,k]
                v_ym = flow.data['vel.y'][i,j,k-1]
                v_zm = flow.data['vel.z'][i,j-1,k]
	        vort_x = 0.5 * ( (v_zp-v_z)/(yp-y) + (v_z-v_zm)/(y-ym) - 2.0*(v_y-v_ym)/(z-zm) )
	    else:
                yp = flow.data['pos.y'][i,j+1,k]
                zp = flow.data['pos.z'][i,j,k+1]
                ym = flow.data['pos.y'][i,j-1,k]
                zm = flow.data['pos.z'][i,j,k-1]
                v_yp = flow.data['vel.y'][i,j,k+1]
                v_zp = flow.data['vel.z'][i,j+1,k]
                v_ym = flow.data['vel.y'][i,j,k-1]
                v_zm = flow.data['vel.z'][i,j-1,k]
	        vort_x = 0.5 * ( (v_zp-v_z)/(yp-y) + (v_z-v_zm)/(y-ym) - (v_yp-v_y)/(zp-z) - (v_y-v_ym)/(z-zm) )
	    Circ += math.fabs(vort_x) * face_area
            vel_abs = Vector(v_x, v_y, v_z)
	    M = math.sqrt(v_x*v_x + v_y*v_y + v_z*v_z)/a
	    gam = 1.4
	    p0 = pressure*math.pow(1.0+0.5*(gam-1.0)*M*M,gam/(gam-1.0))
            df = face_area * pressure * face_normal
            dm_flux = rho * dot(vel_abs, face_normal) * face_area
            m_flux += dm_flux
	    f0 = flow.data['massf[0]-H2'][i,j,k]
            f1 = flow.data['massf[1]-air'][i,j,k]
            mf_flux += f0*dm_flux
            if f0 < 0.02961*f1:
                mfR_flux += f0*dm_flux
            else:
                mfR_flux += 0.02961*f1*dm_flux
            if f0 > f0lim:
                if z > z_pen:
		   z_pen = z
                   zp = flow.data['pos.z'][i,j,k+1]
                   f0p = flow.data['massf[0]-H2'][i,j,k+1]
		   if f0p < f0lim:
		      z_pen = (f0lim-f0)/(f0p-f0)*(zp-z)+z
            prec += p0*dm_flux
            H_flux += dm_flux * (e + pressure/rho + 0.5*math.sqrt(v_x*v_x + v_y*v_y + v_z*v_z))
    return m_flux, H_flux, mf_flux, mfR_flux, prec, Circ, z_pen

fp = open(outputFileName, "w")
fp.write("# 1:x(m) 2:mass_flux(kg/s) 3:total_enthalpy_flux(W) ")
fp.write("4:fuel_flux(kg/s) 5:reactable_fuel_flux(kg/s) 6:pressure recovery(Pa.kg/s) ")
fp.write("7:z_pen(m) 8: Abs_Streamwise_Circulation\n")

for ib in range(nib):
    for i in range(flow[ib][jb].ni):
        m_flux = 0.0; H_flux = 0.0
        mf_flux = 0.0; mfR_flux = 0.0; prec = 0.0; z_pen = 0.0; Circ=0.0
        x = flow[ib][jb].data['pos.x'][i][0][0] # Assuming planar slices
        print x,
        sys.stdout.flush()
        for jb in range(njb):
            m, H, mf, mfR, prec_njb, Circ_njb, zp = \
                process_one_slice(flow[ib][jb], grd[ib][jb], i)
            m_flux += m; H_flux += H
            mf_flux += mf; mfR_flux += mfR; prec += prec_njb; Circ += Circ_njb
	    if zp>z_pen:
	       z_pen = zp
        fp.write("%g %g %g %g %g %g %g %g\n" % 
                 (x, m_flux, H_flux, mf_flux, mfR_flux, prec,
                  z_pen, Circ))

fp.close()
print "Done."


