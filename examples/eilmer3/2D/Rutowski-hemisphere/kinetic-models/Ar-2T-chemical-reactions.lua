-- Ar-2T-chemical-reactions.lua
--
-- Original reaction rates from:
--
-- Hoffert, M.I. and Lien, H. (1967)
-- Quasi-one-dimensional, nonequilibrium gas dynamics of partially 
-- ionized two-temperature Argon
-- Physics of Fluids, Volume 10 Number 8 pp 1769-1777 Aug. 1967
--
-- New cross-sections from:
--
-- Glass, I.I and Liu, W.S. (1978)
-- Effects of hydrogen impurities on shock structure and stability in 
-- ionizing monatomic gases. Part 1. Argon
-- Journal of Fluid Mechanics, Vol. 84 Part 1 pp 55-77 1978
--
-- The presented rates are in the form:
--     k = A ( Ta / T + 2 ) exp( - Ta / T )
-- and have therefore been curve-fitted to the Generalized Arrhenius 
-- form for numerical implementation.  The Argon impact reaction was 
-- curve fitted in the temperature range of [10500,35000] and the 
-- electron impact reaction was curve fitted in the temperature range of
-- [500,20000].  The respective maximum errors were 0.0259% and 0.0078%.
--
-- Author: Daniel F. Potter
-- Date: 18-Apr-2012
-- Place: DLR, Goettingen, Germany
--
-- History: 
--   18-Apr-2012: - Initial implementation
--   01-Oct-2013: - Updated with better cross-sections from Glass and Liu
--   28-May-2014: - Aesthetic improvements

scheme_t = {
    update = "chemical kinetic ODE MC",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

-- Argon-impact ionization

Q_hoffert = 1.2e-19	-- cm2/eV
Q_glass   = 1.0e-19	-- cm2/eV
f         = Q_glass / Q_hoffert

reaction{
   'Ar + Ar <=> Ar+ + Ar + e-',
   fr={'Arrhenius', A=f*8.996906e06, n=1.004, T_a=129441.6 },
   ec={model='from thermo',iT=-1,species="Ar", mode="translation"}
}

-- Electron-impact ionization

Q_hoffert = 7.0e-18	-- cm2/eV
Q_glass   = 4.9e-18	-- cm2/eV
f         = Q_glass / Q_hoffert

reaction{
   'Ar + e- <=> Ar+ + e- + e-',
   fr={'Park', A=f*9.039202e11, n=0.867, T_a=132482.8, 
       p_name='e_minus', p_mode='translation', s_p=1.0, 
       q_name='NA', q_mode='NA'
   },
   chemistry_energy_coupling={ 
      {species='e_minus', mode='translation', 
       model='electron impact ionization', T_I=181700.0}
   },
   ec={model='from thermo',iT=-1,species="e_minus", mode="translation"}
}

