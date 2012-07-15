-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 2nd excited state of atomic Nitrogen, 2s22p4 1S

O_1S = {}
O_1S.M = {
   value = 15.9994e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermO_1S.inp'
}

-- Nonequilibrium data

O_1S.species_type = "monatomic"
O_1S.h_f = {
   value = 1.557402e+07,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermO_1S.inp'
}
O_1S.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
O_1S.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.      n      E(cm-1)      g     l     L     S    parity 
   -- ===========================================================
   ilev_0   =  { 2,    33792.58,     1,   -1,    0,    0,    2 }
   -- ===========================================================
}

