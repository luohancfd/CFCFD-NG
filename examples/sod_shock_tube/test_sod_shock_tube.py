import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from libfpe import *

def main():
    t = 1.0e-3

    # choose your gas model
    #gasfile = "../../input_files/ideal-air.lua"
    gasfile = "thermally-perfect-air.lua"
    #gasfile = "../../input_files/noble-abel-air.lua"
    #gasfile = "../../input_files/van-der-waals-air.lua"

    set_gas_model(gasfile)
    g = get_gas_model_ptr()

    # downstream gas
    Q1 = Gas_data(g)

    molef = {'O2':0.21, 'N2':0.79}
    set_molef(Q1, g, molef)

    Q1.rho = 1.0
    Q1.p = 1e4
    g.eval_thermo_state_rhop(Q1)

    # upstream gas
    Q4 = Gas_data(g)

    molef = {'O2':0.21, 'N2':0.79}
    set_molef(Q4, g, molef)

    Q4.rho = 1.0
    Q4.p = 1e5
    g.eval_thermo_state_rhop(Q4)
    
    test_sod_shock_tube("sod-numerical.dat", t, Q1, Q4)
    
if __name__ == '__main__':
    main()
