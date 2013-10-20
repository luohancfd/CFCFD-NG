-- Filename: argon-radiators-parade.py
--   Author: Daniel F. Potter
--     Date: 20th of October 2013
--    Usage: ready to use as the input file for creating the radiation model
--           e.g. in python:
--           import radpy
--           rsm = create_radiation_spectral_model( "argon-radiators-spradian.lua" )

spectral_data = {
   spectral_model = 'spradian',
   iT = 0,
   iTr = 0,
   iTv = 0,
   iTe = 0,
   radiators = { 'Ar', 'Ar_plus', 'e_minus' },
   lambda_min = 66.666667,
   lambda_max = 10000.000000,
   spectral_points = 14900,
   spectral_blocks = 1,
   adaptive_spectral_grid = false,
}

Ar = {}
Ar.isp = 0
Ar.type = 'atomic radiator'
Ar.mol_weight = 3.994800e-02

Ar_plus = {}
Ar_plus.isp = 0
Ar_plus.type = 'atomic radiator'
Ar_plus.mol_weight = 3.994745e-02

e_minus = {}
e_minus.isp = 0
e_minus.type = 'atomic radiator'
e_minus.mol_weight = 5.485799e-07