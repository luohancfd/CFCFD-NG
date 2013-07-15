-- TC2M5_reactions.lua
--
-- Arrhenius reaction rates from:
--
-- TC2 – M5: Definition of shock tunnel test-cases for gas radiation prediction
-- in Mars-like atmosphere (2012)
-- 
--
-- Author: Daniel F. Potter
-- Date: 15-July-2012
-- Place: Göttingen, Germany
--
-- History: 
--   15-July-2012: First cut
--   15-July-2013: Thermally perfect gas version

scheme_t = {
    update = "chemical kinetic ODE MC",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

-- 3. Heavy particle-impact dissociation

-- 3.1 C2

reaction{
    'C2 + M <=> C + C + M',
    fr={'Arrhenius', A=1.50e16, n=0.0, T_a=71600.0}
}

-- 3.2 CO2

reaction{
    'CO2 + M <=> CO + O + M',
    fr={'Arrhenius', A=1.40e22, n=-1.50, T_a=63275.0},
    efficiencies={CO2=0.0,CO=0.0,CO_plus=0.0,O2=0.0,N2=0.0,NO=0.0,CN=0.0,C2=0.0,C=1.0,C_plus=1.0,O=1.0,O_plus=1.0,N=1.0,N_plus=1.0,e_minus=0.0}
}

reaction{
    'CO2 + M <=> CO + O + M',
    fr={'Arrhenius', A=6.90e21, n=-1.50, T_a=63275.0},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=1.0,O2=1.0,N2=1.0,NO=1.0,CN=1.0,C2=1.0,C=0.0,C_plus=0.0,O=0.0,O_plus=0.0,N=0.0,N_plus=0.0,e_minus=0.0}
}

-- 3.3 CN

reaction{
    'CN + M <=> C + N + M',
    fr={'Arrhenius', A=2.53e14, n=0.0, T_a=88672.0},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=1.0,O2=1.0,N2=1.0,NO=1.0,CN=1.0,C2=1.0,C=1.0,C_plus=1.0,O=1.0,O_plus=1.0,N=1.0,N_plus=1.0,e_minus=0.0}
}

-- 3.3 CO

reaction{
    'CO + M <=> C + O + M',
    fr={'Arrhenius', A=3.40e20, n=-1.00, T_a=128755.0},
    efficiencies={CO2=0.0,CO=0.0,CO_plus=0.0,O2=0.0,N2=0.0,NO=0.0,CN=0.0,C2=0.0,C=1.0,C_plus=0.0,O=1.0,O_plus=0.0,N=1.0,N_plus=0.0,e_minus=0.0}
}

reaction{
    'CO + M <=> C + O + M',
    fr={'Arrhenius', A=2.30e20, n=-1.00, T_a=128755.0},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=0.0,O2=1.0,N2=1.0,NO=1.0,CN=0.0,C2=0.0,C=0.0,C_plus=0.0,O=0.0,O_plus=0.0,N=0.0,N_plus=0.0,e_minus=0.0}
}

-- 3.3 N2

reaction{
    'N2 + M <=> N + N + M',
    fr={'Arrhenius', A=3.00e22, n=-1.60, T_a=113288.0},
    efficiencies={CO2=0.0,CO=0.0,CO_plus=0.0,O2=0.0,N2=0.0,NO=0.0,CN=0.0,C2=0.0,C=1.0,C_plus=0.0,O=1.0,O_plus=0.0,N=1.0,N_plus=0.0,e_minus=0.0}
}

reaction{
    'N2 + M <=> N + N + M',
    fr={'Arrhenius', A=7.01e21, n=-1.60, T_a=113288.0},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=0.0,O2=1.0,N2=1.0,NO=1.0,CN=0.0,C2=0.0,C=0.0,C_plus=0.0,O=0.0,O_plus=0.0,N=0.0,N_plus=0.0,e_minus=0.0}
}

-- 3.4 O2

reaction{
    'O2 + M <=> O + O + M',
    fr={'Arrhenius', A=1.0e22,  n=-1.5, T_a=59750.0},
    efficiencies={CO2=0.0,CO=0.0,CO_plus=0.0,O2=0.0,N2=0.0,NO=0.0,CN=0.0,C2=0.0,C=1.0,C_plus=0.0,O=1.0,O_plus=0.0,N=1.0,N_plus=0.0,e_minus=0.0}
}

reaction{
    'O2 + M <=> O + O + M',
    fr={'Arrhenius', A=2.0e21,  n=-1.5, T_a=59750.0},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=0.0,O2=1.0,N2=1.0,NO=1.0,CN=0.0,C2=0.0,C=0.0,C_plus=0.0,O=0.0,O_plus=0.0,N=0.0,N_plus=0.0,e_minus=0.0}
}

-- 3.5 NO

reaction{
    'NO + M <=> N + O + M',
    fr={'Arrhenius', A=9.63e14,  n=0.0, T_a=75210.0},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=0.0,O2=0.0,N2=0.0,NO=1.0,CN=0.0,C2=0.0,C=1.0,C_plus=0.0,O=1.0,O_plus=0.0,N=1.0,N_plus=0.0,e_minus=0.0}
}

reaction{
    'NO + M <=> N + O + M',
    fr={'Arrhenius', A=1.45e15, n=0.0, T_a=75210.0},
    efficiencies={CO2=0.0,CO=0.0,CO_plus=0.0,O2=1.0,N2=1.0,NO=0.0,CN=0.0,C2=0.0,C=0.0,C_plus=0.0,O=0.0,O_plus=0.0,N=0.0,N_plus=0.0,e_minus=0.0}
}

-- 4. Electron-impact dissociation

-- not in mine, except N2

-- 4.1 CN

reaction{
    'CN + e- <=> C + N + e-',
    fr={'Arrhenius', A=3.05e13, n=0.432, T_a=88966.0}
}

-- 4.2 CO

reaction{
    'CO + e- <=> C + O + e-',
    fr={'Arrhenius', A=8.65e12, n=0.367, T_a=129271.0}
}
-- 4.3 N2

reaction{
   'N2 + e- <=> N + N + e-',
   fr={'Arrhenius', A=2.47e-9, n=6.16, T_a=113263.0}
}

-- 4.4 O2

reaction{
    'O2 + e- <=> O + O + e-',
    fr={'Arrhenius', A=3.47e2,  n=3.52, T_a=59370.0}
}

-- 4.5 NO

reaction{
    'NO + e- <=> N + O + e-',
    fr={'Arrhenius', A=1.05e-2,  n=4.52, T_a=75390.0}
}

-- 5. Neutral exchange

reaction{
   'N2 + O <=> NO + N',
   fr={'Arrhenius', A=5.69e12,  n=0.42, T_a=42938.0}
}

-- not in mine

reaction{
   'O2 + N <=> NO + O',
   fr={'Arrhenius', A=2.49e12,  n=1.179, T_a=4005.0}
}

reaction{
   'CO + C <=> C2 + O',
   fr={'Arrhenius', A=2.0e17,  n=-1.0, T_a=58000.0}
}

reaction{
   'CO + O <=> O2 + C',
   fr={'Arrhenius', A=3.9e13,  n=-0.18, T_a=69200.0}
}

reaction{
   'CO + N <=> CN + O',
   fr={'Arrhenius', A=1.0e14,  n=0.0, T_a=38600.0}
}

reaction{
   'N2 + C <=> CN + N',
   fr={'Arrhenius', A=5.24e13,  n=0.0, T_a=22600.0}
}

reaction{
   'CN + O <=> NO + C',
   fr={'Arrhenius', A=1.6e13,  n=0.1, T_a=14600.0}
}

reaction{
   'CN + C <=> C2 + N',
   fr={'Arrhenius', A=5.0e13,  n=0.0, T_a=13000.0}
}

reaction{
   'CO + CO <=> C + CO2',
   fr={'Arrhenius', A=2.3e9,  n=0.5, T_a=65710.0}
}

reaction{
   'CO2 + O <=> O2 + CO',
   fr={'Arrhenius', A=2.1e13,  n=0.0, T_a=27800.0}
}

reaction{
   'C2 + N2 <=> CN + CN',
   fr={'Arrhenius', A=1.5e13,  n=0.0, T_a=21000.0}
}

reaction{
   'CO + N <=> NO + C',
   fr={'Arrhenius', A=2.9e11,  n=0.5, T_a=53630.0}
}

reaction{
   'NO + CO <=> CO2 + N',
   fr={'Arrhenius', A=4.6e8,  n=0.5, T_a=12070.0}
}

-- not in mine

reaction{
   'NO + NO <=> N2 + O2',
   fr={'Arrhenius', A=3.07e11,  n=0.0, T_a=33660.0}
}

-- 8. Associative ionisation

-- check that the coupling model works when the reaction is written backwards

reaction{
   'CO+ + e- <=> C + O',
   fr={'Arrhenius', A=1.86e18, n=-0.48, T_a=0.0}
}

-- 9. Electron impact ionisation 

reaction{
   'C + e- <=> C+ + e- + e-',
   fr={'Arrhenius', A=5.85e27, n=-2.074, T_a=127510.0}
}

reaction{
   'N + e- <=> N+ + e- + e-',
   fr={'Arrhenius', A=1.93e31, n=-2.856, T_a=168970.0}
}

-- NOTE: the A parameter in the following reaction has been corrected from 8.25e37 to 8.25e27, based on comparison with the park rate

reaction{
   'O + e- <=> O+ + e- + e-',
   fr={'Arrhenius', A=8.25e27, n=-2.237, T_a=157840.0}
}

-- 10. Charge exchange reactions

reaction{
   'CO + C+ <=> CO+ + C',
   fr={'Arrhenius', A=1.0e13, n=0.0, T_a=31400.0}
}

-- not in mine

reaction{
   'CO+ + O <=> CO + O+',
   fr={'Arrhenius', A=8.43e13, n=0.0, T_a=0.0}
}

-- not in mine

reaction{
   'N+ + O <=> N + O+',
   fr={'Arrhenius', A=6.02e11, n=0.0, T_a=0.0}
}

-- not in mine

reaction{
   'O+ + N <=> N+ + O',
   fr={'Arrhenius', A=7.83e13, n=0.0, T_a=0.0}
}

-- 11. Ion exchange reactions

reaction{
   'N+ + O2 <=> O+ + NO',
   fr={'Arrhenius', A=1.69e13, n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N+ + NO <=> N2 + O+',
   fr={'Arrhenius', A=6.02e11, n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

