-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 1st excited state of atomic Nitrogen, 2s22p4 1D

O_1D = {}
O_1D.M = {
   value = 15.9994e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermO_1D.inp'
}

-- Nonequilibrium data

O_1D.species_type = "monatomic"
O_1D.h_f = {
   value = 1.557402e+07,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermO_1D.inp'
}
O_1D.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
O_1D.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.      n      E(cm-1)      g     l     L     S    parity 
   -- ===========================================================
   ilev_0   =  { 2,    15867.86,     5,   -1,    2,    0,    2 }
   -- ===========================================================
}

