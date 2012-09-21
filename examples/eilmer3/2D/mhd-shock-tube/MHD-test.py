## \file MHD-test.py
## \brief Simple setup file for testing MHD implementation in Eilmer3.
## \author DB, 03-Sep-2012

from math import cos, sin, sqrt

job_title = "Simple MHD shock tube test"
print job_title

# We can set individual attributes of the global data object.
gdata.dimensions = 2
gdata.title = job_title
gdata.mhd_flag = 1
gdata.axisymmetric_flag = 0
gdata.stringent_cfl = 1  # to match the old mb_cns behaviour

# Accept defaults for He
select_gas_model(model='ideal gas', species=['He'])

rho0 = 0.1786
p0 = 101.325e3
B0 = sqrt(p0)
R = 2077.0

L_inf = 1.0
u_inf = sqrt(p0/rho0)
t_inf = L_inf/u_inf

rhoL = 3
pL = 3
BxL = 1.5
ByL = 1.0
BzL = 0.0

left = FlowCondition(p=pL*p0, T=(pL*p0)/(rhoL*rho0*R), Bx=BxL*B0, By=ByL*B0, Bz=BzL*B0)

rhoR = 1
pR = 1
BxR = 1.5
ByR = cos(1.5)
BzR = sin(1.5)

right = FlowCondition(p=pR*p0, T=(pR*p0)/(rhoR*rho0*R), Bx=BxR*B0, By=ByR*B0, Bz=BzR*B0)

start_x = -1.0
mid_x = 0.0
stop_x = 1.0

a = Node(mid_x, 0.0, label="a"); b = Node(mid_x, 0.1, label="b")
c = Node(start_x, 0.1, label="c"); d = Node(start_x, 0.0,  label="d")
e = Node(stop_x, 0.0, label="e"); f = Node(stop_x, 0.1, label="f")
south0 = Line(d, a); south1 = Line(a, e) # lower boundary along x-axis
north0 = Line(c, b); north1 = Line(b, f) # upper boundary
# left-end, diaphragm, right-end
west0 = Line(d, c); east0west1 = Line(a, b); east1 = Line(e, f) 

# Define the blocks, boundary conditions and set the discretisation.
blk_0 = Block2D(make_patch(north0, east0west1, south0, west0), 
                nni=800, nnj=2,
                fill_condition=left, label="left")
blk_1 = Block2D(make_patch(north1, east1, south1, east0west1), 
                nni=800, nnj=2,
                fill_condition=right, label="right")
                
identify_block_connections()

blk_0.set_BC(WEST, EXTRAPOLATE_OUT)
blk_1.set_BC(EAST, EXTRAPOLATE_OUT)

# periodic boundary conditions in the y direction
connect_blocks_2D(blk_0, NORTH, blk_0, SOUTH)
connect_blocks_2D(blk_1, NORTH, blk_1, SOUTH)


# Some simulation parameters
gdata.flux_calc = HLLE
gdata.max_time = 0.4*t_inf
gdata.max_step = 8000
gdata.dt = 1.0e-9
#gdata.dt_plot = 1e-4


sketch.xaxis(start_x, stop_x, 0.5, -0.05)
sketch.yaxis(0.0, 0.1, 0.1, -0.05)
sketch.window(start_x, 0.0, stop_x, 0.1, 0.02, 0.02, 0.17, 0.035)


# EXACT SOLUTION
import numpy
t=0.4
xexact=numpy.array([-2.500000, -1.474922, -0.990247, -0.631585, -0.631585, 
                 -0.521395, -0.445268, 0.402052, 0.402052, 1.279598, 
                 1.279598,  1.568067, 1.568067, 2.072332, 2.072332, 
                 2.500000])*t
rexact=numpy.array([ 3.000000, 3.000000, 2.340949, 2.340949, 2.340949, 2.340949, 
        2.200167, 2.200167, 1.408739, 1.408739, 1.054703, 1.054703, 1.054703, 
        1.054703, 1.000000, 1.000000])*rho0
Byexact=numpy.array([1.000000, 1.000000, 0.642777, 0.642777, 0.344252, 0.344252, 
         0.413199, 0.413199, 0.413199, 0.413199, 0.601050, 0.601050, 
         0.079386, 0.079386, 0.070737, 0.070737])*B0
Bzexact=numpy.array([0.000000, 0.000000, 0.000000, 0.000000, 0.542820, 0.542820, 
         0.651535, 0.651535, 0.651535, 0.651535, 0.947741, 0.947741, 
         1.119452, 1.119452, 0.997495, 0.997495])*B0

f = open("exact-sol.dat","w")
f.write("# Variables: 1:x 2:density 3:B.y 4:B.z\n")
for i in range(len(rexact)):
    f.write("%1.8f %1.8f %1.8f %1.8f\n"%(xexact[i],rexact[i],Byexact[i],Bzexact[i]))

f.close()
