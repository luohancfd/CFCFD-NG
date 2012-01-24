# file: lehr1.py
#
# Spherical nose of Lehr's projectile in detonable gas -- continued.
#
# PJ, 27-Feb-2010
# Adapted bits from sphere-heat-transfer and Rowan's mbcsn2/lehr_sphere.
# This is the continuation of the simulation on a better fitted grid.

gdata.title = "Lehr experiment M=4.79"
R = 7.5e-3  # Nose radius, metres

p_inf = 320.0/760.0*101325.0  # Pascals
u_inf = 1931  # m/s
T_inf = 292 # degrees K
p_init = p_inf / 5

select_gas_model(model='thermally perfect gas', 
                 species=['O2','N2','H2','O','H','H2O','OH','HO2'])
# species index             0    1    2   3   4     5    6     7 
set_reaction_scheme("Evans_Schexnayder.lua",reacting_flag=1)

# Calculation: convert mole fractions to mass fractions.
MW_O2 = 3.19988000e-02  # kg/mole
MW_N2 = 2.80134800e-02
MW_H2 = 2.01588000e-03
# moles for a stoichiometric mix
m_O2 = 1.0; m_N2 = 3.76; m_H2 = 2.0
# mole fractions
mole_tot = m_O2 + m_N2 + m_H2
X_O2 = m_O2 / mole_tot
X_N2 = m_N2 / mole_tot
X_H2 = m_H2 / mole_tot
MW_mix = X_O2 * MW_O2 + X_N2 * MW_N2 + X_H2 * MW_H2
# mass fractions
mf = {'O2':X_O2*(MW_O2/MW_mix), 
      'N2':X_N2*(MW_N2/MW_mix), 
      'H2':X_H2*(MW_H2/MW_mix)}
print "mass fractions=", mf

inflow = FlowCondition(p=p_inf, u=u_inf, T=T_inf, massf=mf)
initial = ExistingSolution('lehr', '.', 24, 9999)

# Job-control information
t_final = 50.0e-6
ni = 200; nj = 300
gdata.axisymmetric_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = t_final
gdata.max_step = 800000
gdata.dt = 1.0e-9
gdata.cfl = 0.40
gdata.dt_plot = 1.0e-6
gdata.dt_history = 0.01e-6  # want to capture MHz frequency
 
# Begin geometry details for a single region around a spherical nose.
# The node coordinates are scaled with the body radius.
a = Node(0.0, 0.0, label="a")
b = Node(-1.0*R, 0.0, label="b")
c = Node(0.0, R, label="c")
# The inflow boundary is a Bezier curve.
d = [Node(-1.3*R,0), Node(-1.3*R,0.7*R), Node(-0.87*R,1.4*R), Node(0,2.1*R)] 
# order of boundaries: N, E, S, W
flow_domain0 = make_patch(Line(d[-1],c), Arc(b,c,a), Line(d[0],b), Bezier(d))
boundary_conditions0 = [ExtrapolateOutBC(), SlipWallBC(),
                        SlipWallBC(), SupInBC(inflow)]
blk = SuperBlock2D(psurf=flow_domain0, fill_condition=initial,
                   nni=ni, nnj=nj, nbi=4, nbj=6,
                   bc_list=boundary_conditions0,
                   label="blk")
HistoryLocation(-R,0.0)
HistoryLocation(-R,0.001)

sketch.xaxis(-15.0e-3, 5.0e-3, 5.0e-3, -0.002)
sketch.yaxis(0.0, 20.0e-3, 5.0e-3, 0.0)
sketch.window(-1.5*R, 0.0, 1.5*R, 3.0*R, 0.05, 0.05, 0.15, 0.15)

