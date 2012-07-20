-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 2nd excited state of atomic Nitrogen, 2s22p3 2P¡

N_2Po = {}
N_2Po.M = {
   value = 14.0067e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}

-- Nonequilibrium data

N_2Po.species_type = "monatomic"
N_2Po.h_f = {
   value = 3.374671e+07,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N_2Po.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N_2Po.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.     n      E(cm-1)      g     l     L     S    parity 
   -- ===========================================================
   ilev_0   =  { 2,    28839.18,     6,   -1,    1,    1,    1 }
   -- ===========================================================
