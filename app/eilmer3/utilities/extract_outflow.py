#!/usr/bin/env python
# extract_outflow.py
#
# Custom postprocessor to extract the outflow plane of a given simulation (jobname input)
# map it onto a new grid in a mass flux weighted manner and write out a table that can
# be input into a user defined boundary condition.  Fully conservative mapping is not
# possible due to the EOS being unknown.
#
# Must be executed in the source simulation directory with the following arguments:
#
# extract_outflow.py SourceJobName DestinationJobName DestinationDirectory outputFileName
#
# Example of usage:
#
# extract_outflow.py 5_enriched_inject 5_enriched_inject ../../ArnoldMafu/EP05/ inflow.lua
#
# PJ, 21-May-2010 & VW, 1-August-2013

import sys, os, gzip, math
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_grid import StructuredGrid
from e3_flow import StructuredGridFlow
from libprep3 import Vector, quad_centroid, quad_normal, quad_area, cross, dot

if len(sys.argv) != 5:
    print "Usage (in source simulation directory): extract_outflow.py SourceJobName DestinationJobName DestinationDirectory outputFileName"
    sys.exit()
jobName = sys.argv[1]
destJobName = sys.argv[2]
dest_dir = sys.argv[3]
outputFileName = sys.argv[4]
tindx = '9999'

output_format = 1	# 1 = lua table, 2 = tecplot file

# SOURCE SIMULATION
# We'll access the blocks at the downstream end of the grid using [j][k] indexing, 
njb = 4
nkb = 1
bindx_list = [[14],[21],[22],[23]]
flow = []
grd = []
jmax=0
kmax=0
for jb in range(njb):
    flow.append([])
    grd.append([])
    for kb in range(nkb):
        bindx = '%04d' % bindx_list[jb][kb]
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
	if jb==0:
		kmax=kmax+g.nk
	if kb==0:
		jmax=jmax+g.nj

# Convert outflow plan to a single patch grid and construct vectors of cell edge locations
ymin = []
ymax = []
zmin = []
zmax = []
y = eval(`[[0]*kmax]*jmax`)
z = eval(`[[0]*kmax]*jmax`)
A = eval(`[[0]*kmax]*jmax`)
rho = eval(`[[0]*kmax]*jmax`)
p = eval(`[[0]*kmax]*jmax`)
T = eval(`[[0]*kmax]*jmax`)
v_x = eval(`[[0]*kmax]*jmax`)
v_y = eval(`[[0]*kmax]*jmax`)
v_z = eval(`[[0]*kmax]*jmax`)
tke = eval(`[[0]*kmax]*jmax`)
omega = eval(`[[0]*kmax]*jmax`)
f0 = eval(`[[0]*kmax]*jmax`)
f1 = eval(`[[0]*kmax]*jmax`)
f2 = eval(`[[0]*kmax]*jmax`)
gamma = eval(`[[0]*kmax]*jmax`)
jglobal = -1
for jb in range(njb):
    i=grd[jb][0].ni-2
    for j in range(grd[jb][0].nj-1):
		kglobal = -1
		jglobal = jglobal+1
        	vtx = grd[jb][0].get_vertex_list_for_cell(i, j, 0)
		# The EAST cell face has 1, 2, 6, 5 as corners (unit normal out).
		#print "jglobal=",jglobal,"vtx[1].y=",vtx[1].y
		ymin.append(vtx[1].y)
		ymax.append(vtx[2].y)
		for kb in range(nkb):
			for k in range(grd[jb][kb].nk-1):
				kglobal = kglobal+1
            			vtx = grd[jb][kb].get_vertex_list_for_cell(i, j, k)
				if jb==0:
					if j==0:
						zmin.append(vtx[1].z)
						zmax.append(vtx[5].z)
				y[jglobal][kglobal] = flow[jb][kb].data['pos.y'][i,j,k]
				z[jglobal][kglobal] = flow[jb][kb].data['pos.z'][i,j,k]
				A[jglobal][kglobal] = quad_area(vtx[1], vtx[2], vtx[6], vtx[5])
				rho[jglobal][kglobal] = flow[jb][kb].data['rho'][i,j,k] 
				p[jglobal][kglobal] = flow[jb][kb].data['p'][i,j,k] 
				T[jglobal][kglobal] = flow[jb][kb].data['T[0]'][i,j,k]
				v_x[jglobal][kglobal] = flow[jb][kb].data['vel.x'][i,j,k]
				v_y[jglobal][kglobal] = flow[jb][kb].data['vel.y'][i,j,k]
				v_z[jglobal][kglobal] = flow[jb][kb].data['vel.z'][i,j,k]
				tke[jglobal][kglobal] = flow[jb][kb].data['tke'][i,j,k]
				omega[jglobal][kglobal] = flow[jb][kb].data['omega'][i,j,k]
				f0[jglobal][kglobal] = flow[jb][kb].data['massf[0]-H2'][i,j,k]
				f1[jglobal][kglobal] = flow[jb][kb].data['massf[1]-air'][i,j,k]
				f2[jglobal][kglobal] = flow[jb][kb].data['massf[2]-O2'][i,j,k]
				gamma[jglobal][kglobal] = flow[jb][kb].data['k[0]'][i,j,k]
				# Note: after this loop is complete jglobal and kglobal store max indices

# DESTINATION GRID
# Assume this is a single block at present 
dnjb = 1
dnkb = 1
dbindx_list = [[0]]
dgrd = []
for jb in range(dnjb):
    dgrd.append([])
    for kb in range(dnkb):
        bindx = '%04d' % dbindx_list[jb][kb]
        print "destination bindx=", bindx
        # The grid is always at tindx 0.
        fileName = dest_dir+'grid/t0000/'+ destJobName + '.grid.b'+bindx+'.t0000.gz'
        fp = gzip.open(fileName, "r")
        dg = StructuredGrid()
        dgrd[-1].append(dg)
        dg.read(fp)
        fp.close()
        print "destination grid data: ni,nj,nk=", dg.ni, dg.nj, dg.nk

# Construct vectors of cell edge locations (assuming rectilinear grids):
dymin = []
dymax = []
dzmin = []
dzmax = []
for jb in range(dnjb):
    for j in range(dgrd[jb][0].nj-1):
        vtx = dgrd[jb][0].get_vertex_list_for_cell(0, j, 0)
	# The EAST cell face has 1, 2, 6, 5 as corners (unit normal out).
	dymin.append(vtx[1].y)
	dymax.append(vtx[2].y)
for kb in range(dnkb):
    for k in range(dgrd[0][kb].nk-1):
        vtx = dgrd[0][kb].get_vertex_list_for_cell(0, 0, k)
	# The EAST cell face has 1, 2, 6, 5 as corners (unit normal out).
	dzmin.append(vtx[1].z)
	dzmax.append(vtx[5].z)		
	
# Determine which source cells overlap with destination cells and compute weights
# Storing both upper and lower limits separate allows boundaries and large source cells to be handled easily
jl = []   		# lowest j index source cell overlapping destination cell dj
ju = []   		# highest j index source cell overlapping destination cell dj
jlweight = []   # fraction (in y) of cell jl overlapping cell dj
juweight = []	# fraction (in y) of cell ju overlapping cell dj
jl.append(0)
for dj in range(dgrd[0][0].nj-1):
	foundedge=0
	j=jl[dj]
	jlweight.append((ymax[j]-dymin[dj])/(ymax[j]-ymin[j]))
	while foundedge==0:
	    	if ymax[j]>=dymax[dj]:
		    	ju.append(j)
			juweight.append((dymax[dj]-ymin[j])/(ymax[j]-ymin[j]))
			# note, if jl=ju, must take into account both weights when mapping
			jl.append(j)
			foundedge=1
		else:
			j=j+1
		if j>jglobal:
			print "Failed to find edge in j direction."
			sys.exit()
			
		

kl = []   		# lowest k index source cell overlapping destination cell dk
ku = []   		# highest k index source cell overlapping destination cell dk
klweight = []   # fraction (in z) of cell kl overlapping cell dk
kuweight = []	# fraction (in z) of cell ku overlapping cell dk
kl.append(0)	
for dk in range(dgrd[0][0].nk-1):
	foundedge=0
	k=kl[dk]
	klweight.append((zmax[k]-dzmin[dk])/(zmax[k]-zmin[k]))
	while foundedge==0:
	    if zmax[k]>=dzmax[dk]:
		ku.append(k)
		kuweight.append((dzmax[dk]-zmin[k])/(zmax[k]-zmin[k]))
		# note, if jl=ju, must take products of weights when mapping
		kl.append(k)
		foundedge=1
	    else:
		k=k+1
		if k>kglobal:
			print "Failed to find edge in k direction."
			sys.exit()

		
# Carry out the mapping and write out the results:
# Note - mass flux weighted averages are used as a fully conservative computation
# of the mapped p and T is not possible due to unknown EOS
print "Mapping and writing data"
fp = open(outputFileName, "w")
jsize=dgrd[0][0].nj-1
ksize=dgrd[0][0].nk-1

if output_format == 1:
	fp.write("prof={}\n")
else:
	fp.write("Variables=,y,z,p,u,v,w,T,tke,omega,f0,f1,f2\n")
	fp.write("Zone i=%d j=%d k=%d\n" % (1,jsize,ksize))
	
for dj in range(dgrd[0][0].nj-1):
	for dk in range(dgrd[0][0].nk-1):
            	vtx = dgrd[0][0].get_vertex_list_for_cell(0, dj, dk)
            	face_centroid = quad_centroid(vtx[1], vtx[2], vtx[6], vtx[5])
		dy = face_centroid.y
		dz = face_centroid.z
		dp = 0
		du = 0
		dv = 0
		dw = 0
		dT = 0
		dtke = 0
		domega = 0
		df0 = 0
		df1 = 0
		df2 = 0
		mtot = 0
		for j in range(jl[dj], ju[dj]):
			if jl[dj]==ju[dj]:
				jweight = jlweight[dj] + juweight[dj] - 1
			elif j==jl[dj]:
				jweight = jlweight[dj]
			elif j==ju[dj]:
				jweight = juweight[dj]
			else:
				jweight = 1
			for k in range(kl[dk], ku[dk]):
				if kl[dk]==ku[dk]:
					kweight = klweight[dk] + kuweight[dk] - 1
				elif k==kl[dk]:
					kweight = klweight[dk]
				elif k==ku[dk]:
					kweight = kuweight[dk]
				else:
					kweight = 1
				dm = jweight*kweight*rho[j][k]*v_x[j][k]*A[j][k]
				dp =  	  dp     + dm*p[j][k]
				du =      du     + dm*v_x[j][k]
				dv =      dv     + dm*v_y[j][k]
				dw =      dw     + dm*v_z[j][k]
				dT =      dT     + dm*T[j][k]
				dtke =    dtke   + dm*tke[j][k]
				domega =  domega + dm*omega[j][k]
				df0 =     df0    + dm*f0[j][k]
				df1 =     df1    + dm*f1[j][k]
				df2 =     df2    + dm*f2[j][k]
				mtot = mtot+dm
		dp = dp/mtot
		du = du/mtot
		dv = dv/mtot
		dw = dw/mtot
		dT = dT/mtot
		dtke = dtke/mtot
		domega = domega/mtot
		df0 = df0/mtot
		df1 = df1/mtot
		df2 = df2/mtot
		if output_format == 1:        	
			fp.write("prof[%d]={y=%g, z=%g, p=%g, u=%g, v=%g, w=%g, T=%g, tke=%g, omega=%g, f0=%g, f1=%g, f2=%g}\n" % (dj*jsize+dk+1,dy,dz,dp,du,dv,dw,dT,dtke,domega,df0,df1,df2))
		else:
			fp.write("%g %g %g %g %g %g %g %g %g %g %g %g\n" % (dy,dz,dp,du,dv,dw,dT,dtke,domega,df0,df1,df2))

fp.close()
print "Done."


