# Mass-flux boundary-condition demonstration -- SubsonicInBC version.
# PJ, 04-Feb-2015

gdata.title = "One-dimensional tube with subsonic air."
gdata.dimensions = 2
gdata.axisymmetric_flag = True

select_gas_model(model='ideal gas', species=['air',])
initial_air = FlowCondition(p=100.0e3, u=0.0, v=0.0, T=300.0)

p01 = Vector(0.0,0.005); p11 = Vector(0.1, 0.005)
p00 = Vector(0.0,  0.0); p10 = Vector(0.1,   0.0)
blk = Block2D(CoonsPatch(p00,p10,p11,p01),
              nni=10, nnj=2,
              bc_list=[SlipWallBC(), 
                       ExtrapolateOutBC(),
                       SlipWallBC(),
                       SubsonicInBC(initial_air, mass_flux=12.0, assume_ideal=True)],
              fill_condition=initial_air)

gdata.flux_calc = AUSMDV
gdata.max_time = 500.0e-3
gdata.max_step = 200000
gdata.dt = 1.0e-6


