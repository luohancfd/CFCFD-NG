-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- Ground state of atomic Nitrogen, 2s22p3 4S¡

N_4So = {}
N_4So.M = {
   value = 14.0067e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}

-- Nonequilibrium data

N_4So.species_type = "monatomic"
N_4So.h_f = {
   value = 3.374671e+07,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N_4So.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N_4So.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.     n      E(cm-1)      g     l     L     S    parity 
   -- ===========================================================
   ilev_0   =  { 2,        0.00,     4,   -1,    0,    2,    1 }
   -- ===========================================================
