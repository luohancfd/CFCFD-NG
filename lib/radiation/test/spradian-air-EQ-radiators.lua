-- File: spradian-air-EQ-radiators.lua

spectral_data = {
   spectral_model = 'spradian',
   iT = 0,
   iTr = 0,
   iTv = 0,
   iTe = 0,
   radiators = { 'N', 'N_plus', 'O', 'O_plus', 'e_minus', },
   lambda_min = 50.000000,
   lambda_max = 1000.000000,
   spectral_points = 99500,
   spectral_blocks = 1,
   adaptive_spectral_grid = false
}

N2 = {}
N2.isp = 0
N2.type = 'diatomic radiator'
N2.mol_weight = 2.801348e-02

N2_plus = {}
N2_plus.isp = 1
N2_plus.type = 'diatomic radiator'
N2_plus.mol_weight = 2.801293e-02

NO = {}
NO.isp = 2
NO.type = 'diatomic radiator'
NO.mol_weight = 3.000610e-02

O2 = {}
O2.isp = 4
O2.type = 'diatomic radiator'
O2.mol_weight = 3.199880e-02

O2_plus = {}
O2_plus.isp = 5
O2_plus.type = 'diatomic radiator'
O2_plus.mol_weight = 31.9982514e-3

N = {}
N.isp = 6
N.type = 'atomic radiator'
N.mol_weight = 1.400674e-02

N_plus = {}
N_plus.isp = 7
N_plus.type = 'atomic radiator'
N_plus.mol_weight = 1.400619e-02

O = {}
O.isp = 8
O.type = 'atomic radiator'
O.mol_weight = 1.599940e-02

O_plus = {}
O_plus.isp = 9
O_plus.type = 'atomic radiator'
O_plus.mol_weight = 1.599885e-02

e_minus = {}
e_minus.isp = 10
e_minus.type = 'electron radiator'
e_minus.mol_weight = 5.485799e-07
