-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 1st excited state of atomic Carbon, 2s22p2 1D

C_1D = {}
C_1D.M = { 
   value = 12.0107000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}

-- Nonequilibrium data

C_1D.species_type = "monatomic"
C_1D.h_f = {
   value = 59670127.47,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C_1D.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C_1D.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.      n       E(cm-1)     g     l     L     S     parity 
   -- ===========================================================
   ilev_0   =  { 2,    10192.63,     5,   -1,    2,    0,    2 }
   -- ===========================================================
}

