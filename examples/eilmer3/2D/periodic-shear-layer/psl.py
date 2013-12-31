# psl.py
gdata.title = "Periodic shear layer"

select_gas_model(model='ideal gas', species=['He', 'air'])
gdata.diffusion_model = "FicksFirstLaw"
gdata.diffusion_flag = 1

H = 0.010 # layer thickness in metres
L = 0.100 # wavelength in metres

def initial_gas(x, y, z):
    """
    Top and bottom layer of different gases with a basic velocity shear 
    across the interface, plus a streamwise-periodic perturbation 
    that decays away from the interface.
    """
    global H, L
    from math import sin, exp, pi
    p = 100.0e3
    T = 300.0
    U0 = 1000.0
    #
    # The top and bottom streams.
    if y < 0.0:
        massf = {'He':0.1, 'air':0.9}
    else:
        massf = {'air':1.0}
    #
    # The basic velocity shear.
    if y < -H:
        u = -U0
    elif y < H:
        u = y/H * U0
    else:
        u = U0
    #
    # Add perturbation
    V0 = 50.0
    v = V0 * sin(x/L*pi) * exp(-abs(y)/H)
    flow = FlowCondition(p=p, T=T, u=u, v=v, massf=massf, add_to_list=0)
    return flow.to_dict()

#
# Geometry
ymin = -15.0 * H
ymax = 15.0 * H
xmin = -L
xmax = L

a0 = Node(xmin, ymin); a1 = Node(xmin, ymax)
b0 = Node(xmax, ymin); b1 = Node(xmax, ymax)
domain = make_patch(Line(a1,b1), Line(b0,b1), Line(a0,b0), Line(a0,a1))

nnx = 150; nny = 300
nbi = 2; nbj = 2

superblk = SuperBlock2D(psurf=domain, nni=nnx, nnj=nny,
                        bc_list=[SlipWallBC(),]*4,
                        fill_condition=initial_gas,
                        nbi=nbi, nbj=nbj, label="blk")
# Make the domain periodic in the x-direction.
for j in range(nbj):
    connect_blocks_2D(superblk.blks[-1][j], EAST, 
                      superblk.blks[0][j], WEST,
                      check_corner_locations=False)

gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 15.0e-3  # seconds
gdata.max_step = 150000
gdata.dt = 1.0e-9
gdata.dt_plot = gdata.max_time / 50.0

sketch.xaxis(-0.1, 0.1, 0.05, -0.03)
sketch.yaxis(-0.15, 0.15, 0.05, -0.03)
sketch.window(-0.15, -0.15, 0.15, 0.15, 0.05, 0.05, 0.17, 0.17)

