-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 2nd excited state of atomic Carbon, 2s22p2 1S

C_1S = {}
C_1S.M = { 
   value = 12.0107000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}

-- Nonequilibrium data

C_1S.species_type = "monatomic"
C_1S.h_f = {
   value = 59670127.47,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C_1S.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C_1S.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.      n       E(cm-1)     g     l     L     S     parity 
   -- ===========================================================
   ilev_0   =  { 2,    21648.01,     1,   -1,    0,    0,    2 }
   -- ===========================================================
}

