-- File: grey-gas-DT-A.lua
-- NOTE: this file was created by hand, do not delete!

spectral_data = {
   spectral_model = 'photaura',
   radiators = { 'air', },
   lambda_min = 10.000000,
   lambda_max = 3000.000000,
   spectral_points = 501,
   spectral_blocks = 1,
}

transport_data = {
   transport_model = 'monte carlo',
   spectrally_resolved = 1,
   nrays = 512,
   clustering = "none",
   absorption = "partitioned energy"
}

air = {}
air.isp = 0
air.type = 'planck_radiator'
air.E_pop_method = 'NA'
air.mol_weight = 0.028964
air.h_f = 0.0
air.eta_I = 0.0
air.Z = 0
air.iT = 0
air.iTe = 0
air.kappa_const = 1.0


