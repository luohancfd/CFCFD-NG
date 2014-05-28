-- Ar-2T-energy-exchange.lua
--
-- Electron-translation thermal energy exchange for the Ar,Ar+,e- system
-- via the Appleton and Bray (1967) model. Heavy-particle excitation 
-- cross sections have been curve fitted from the data presented in:
--
-- Hoffert, M.I. and Lien, H. (1967)
-- Quasi-one-dimensional, nonequilibrium gas dynamics of partially 
-- ionized two-temperature Argon
-- Physics of Fluids, Volume 10 Number 8 pp 1769-1777 Aug. 1967
--
-- Author: Daniel F. Potter
-- Date: 18-Apr-2012
-- Place: DLR, Goettingen, Germany
--
-- History: 
--   18-Apr-2012: - Initial implementation
--   28-May-2014: - Aesthetic improvements

mechanism{
   'e- ~~ Ar : E-T',
   rt={'Appleton-Bray:TwoRangeNeutral',
        T_switch=10000.0,
	sigma_low_T={  3.9e-21, -5.51e-25, 5.95e-29},
	sigma_high_T={-3.5e-21,  7.75e-25, 0.0}
   }
}

mechanism{
   'e- ~~ Ar+ : E-T',
   rt={'Appleton-Bray:Ion'}
}

