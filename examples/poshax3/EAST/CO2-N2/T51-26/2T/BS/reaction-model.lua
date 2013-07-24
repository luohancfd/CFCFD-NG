-- TC2M5-15sp.lua
--
-- Arrhenius reaction rates from:
--
-- TC2 – M5: Definition of shock tunnel test-cases for gas radiation prediction
-- in Mars-like atmosphere (2012)
-- 
-- Ground state rates have been used where the reaction is not for a bulk-species

--
-- Author: Daniel F. Potter
-- Date: 30-Sep-2012
-- Place: Göttingen, Germany
--
-- History: 
--   30-Sep-2012: Verification
--   11-Oct-2012: Added reactions from Park involving NO+ and O2+

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
    fr={'Park', A=1.50e16, n=0.0, T_a=71600.0, p_name='C2', p_mode='translation', s_p=0.5, q_name='C2', q_mode='vibration'},
    chemistry_energy_coupling={ {species='C2', mode='vibration', model='TreanorMarrone', U=71600.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

-- 3.2 CO2

reaction{
    'CO2 + M <=> CO + O + M',
    fr={'Park', A=1.40e22, n=-1.50, T_a=63275.0, p_name='CO2', p_mode='translation', s_p=0.5, q_name='CO2', q_mode='vibration'},
    efficiencies={CO2=0.0,CO=0.0,CO_plus=0.0,O2=0.0,N2=0.0,NO=0.0,CN=0.0,C2=0.0,
                  C=1.0,C_plus=1.0,O=1.0,O_plus=1.0,N=1.0,N_plus=1.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO2', mode='vibration', model='TreanorMarrone', U=63275.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

reaction{
    'CO2 + M <=> CO + O + M',
    fr={'Park', A=6.90e21, n=-1.50, T_a=63275.0, p_name='CO2', p_mode='translation', s_p=0.5, q_name='CO2', q_mode='vibration'},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=1.0,O2=1.0,N2=1.0,NO=1.0,CN=1.0,C2=1.0,
                  C=0.0,C_plus=0.0,O=0.0,O_plus=0.0,N=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO2', mode='vibration', model='TreanorMarrone', U=63275.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

-- 3.3 CN

reaction{
    'CN + M <=> C + N + M',
    fr={'Park', A=2.53e14, n=0.0, T_a=88672.0, p_name='CN', p_mode='translation', s_p=0.5, q_name='CN', q_mode='vibration'},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=1.0,O2=1.0,N2=1.0,NO=1.0,CN=1.0,C2=1.0,
                  C=1.0,C_plus=1.0,O=1.0,O_plus=1.0,N=1.0,N_plus=1.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CN', mode='vibration', model='TreanorMarrone', U=88672.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

-- 3.3 CO

reaction{
    'CO + M <=> C + O + M',
    fr={'Park', A=3.40e20, n=-1.00, T_a=128755.0, p_name='CO', p_mode='translation', s_p=0.5, q_name='CO', q_mode='vibration'},
    efficiencies={CO2=0.0,CO=0.0,CO_plus=0.0,O2=0.0,N2=0.0,NO=0.0,CN=0.0,C2=0.0,
                  C=1.0,C_plus=0.0,O=1.0,O_plus=0.0,N=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO', mode='vibration', model='TreanorMarrone', U=128755.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

reaction{
    'CO + M <=> C + O + M',
    fr={'Park', A=2.30e20, n=-1.00, T_a=128755.0, p_name='CO', p_mode='translation', s_p=0.5, q_name='CO', q_mode='vibration'},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=0.0,O2=1.0,N2=1.0,NO=1.0,CN=0.0,C2=0.0,
                  C=0.0,C_plus=0.0,O=0.0,O_plus=0.0,N=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO', mode='vibration', model='TreanorMarrone', U=128755.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

-- 3.3 N2

reaction{
    'N2 + M <=> N + N + M',
    fr={'Park', A=3.00e22, n=-1.60, T_a=113288.0, p_name='N2', p_mode='translation', s_p=0.5, q_name='N2', q_mode='vibration'},
    efficiencies={CO2=0.0,CO=0.0,CO_plus=0.0,O2=0.0,N2=0.0,NO=0.0,CN=0.0,C2=0.0,
                  C=1.0,C_plus=0.0,O=1.0,O_plus=0.0,N=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='N2', mode='vibration', model='TreanorMarrone', U=113288.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

reaction{
    'N2 + M <=> N + N + M',
    fr={'Park', A=7.01e21, n=-1.60, T_a=113288.0, p_name='N2', p_mode='translation', s_p=0.5, q_name='N2', q_mode='vibration'},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=0.0,O2=1.0,N2=1.0,NO=1.0,CN=0.0,C2=0.0,
                  C=0.0,C_plus=0.0,O=0.0,O_plus=0.0,N=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='N2', mode='vibration', model='TreanorMarrone', U=113288.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

-- 3.4 O2

reaction{
    'O2 + M <=> O + O + M',
    fr={'Park', A=1.0e22,  n=-1.5, T_a=59750.0, p_name='O2', p_mode='translation', s_p=0.5, q_name='O2', q_mode='vibration'},
    efficiencies={CO2=0.0,CO=0.0,CO_plus=0.0,O2=0.0,N2=0.0,NO=0.0,CN=0.0,C2=0.0,
                  C=1.0,C_plus=0.0,O=1.0,O_plus=0.0,N=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='O2', mode='vibration', model='TreanorMarrone', U=59750.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

reaction{
    'O2 + M <=> O + O + M',
    fr={'Park', A=2.0e21,  n=-1.5, T_a=59750.0, p_name='O2', p_mode='translation', s_p=0.5, q_name='O2', q_mode='vibration'},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=0.0,O2=1.0,N2=1.0,NO=1.0,CN=0.0,C2=0.0,
                  C=0.0,C_plus=0.0,O=0.0,O_plus=0.0,N=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='O2', mode='vibration', model='TreanorMarrone', U=59750.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

-- 3.5 NO

reaction{
    'NO + M <=> N + O + M',
    fr={'Park', A=9.63e14,  n=0.0, T_a=75210.0, p_name='NO', p_mode='translation', s_p=0.5, q_name='NO', q_mode='vibration'},
    efficiencies={CO2=1.0,CO=1.0,CO_plus=0.0,O2=0.0,N2=0.0,NO=1.0,CN=0.0,C2=0.0,
                  C=1.0,C_plus=0.0,O=1.0,O_plus=0.0,N=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='NO', mode='vibration', model='TreanorMarrone', U=75210.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

reaction{
    'NO + M <=> N + O + M',
    fr={'Park', A=1.45e15, n=0.0, T_a=75210.0, p_name='NO', p_mode='vibration', s_p=0.5, q_name='NO', q_mode='translation'},
    efficiencies={CO2=0.0,CO=0.0,CO_plus=0.0,O2=1.0,N2=1.0,NO=0.0,CN=0.0,C2=0.0,
                  C=0.0,C_plus=0.0,O=0.0,O_plus=0.0,N=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='NO', mode='vibration', model='TreanorMarrone', U=75210.0/3.0 } },
    ec={model='from CEA curves',iT=0}
}

-- 4. Electron-impact dissociation

-- 4.1 CN

reaction{
    'CN + e- <=> C + N + e-',
    fr={'Park', A=3.05e13, n=0.432, T_a=88966.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='CN', q_mode='vibration'},
    ec={model='from CEA curves',iT=0}
}

-- 4.2 CO

reaction{
    'CO + e- <=> C + O + e-',
    fr={'Park', A=8.65e12, n=0.367, T_a=129271.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='CO', q_mode='vibration'},
    ec={model='from CEA curves',iT=0}
}

-- 4.3 N2

-- the A parameter for this reaction seems very low, but note that n is very large

reaction{
    'N2 + e- <=> N + N + e-',
    fr={'Park', A=2.47e-9, n=6.16, T_a=113263.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='N2', q_mode='vibration'},
    ec={model='from CEA curves',iT=0}
}

-- 4.4 O2

reaction{
    'O2 + e- <=> O + O + e-',
    fr={'Park', A=3.47e2,  n=3.52, T_a=59370.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='O2', q_mode='vibration'},
    ec={model='from CEA curves',iT=0}
}

-- 4.5 NO

reaction{
    'NO + e- <=> N + O + e-',
    fr={'Park', A=1.05e-2,  n=4.52, T_a=75390.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='NO', q_mode='vibration'},
    ec={model='from CEA curves',iT=0}
}

-- 5. Neutral exchange

reaction{
   'N2 + O <=> NO + N',
   fr={'Arrhenius', A=5.69e12,  n=0.42, T_a=42938.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O2 + N <=> NO + O',
   fr={'Arrhenius', A=2.49e12,  n=1.179, T_a=4005.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + C <=> C2 + O',
   fr={'Arrhenius', A=2.0e17,  n=-1.0, T_a=58000.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + O <=> O2 + C',
   fr={'Arrhenius', A=3.9e13,  n=-0.18, T_a=69200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + N <=> CN + O',
   fr={'Arrhenius', A=1.0e14,  n=0.0, T_a=38600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + C <=> CN + N',
   fr={'Arrhenius', A=5.24e13,  n=0.0, T_a=22600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CN + O <=> NO + C',
   fr={'Arrhenius', A=1.6e13,  n=0.1, T_a=14600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CN + C <=> C2 + N',
   fr={'Arrhenius', A=5.0e13,  n=0.0, T_a=13000.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + CO <=> C + CO2',
   fr={'Arrhenius', A=2.3e9,  n=0.5, T_a=65710.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO2 + O <=> O2 + CO',
   fr={'Arrhenius', A=2.1e13,  n=0.0, T_a=27800.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'C2 + N2 <=> CN + CN',
   fr={'Arrhenius', A=1.5e13,  n=0.0, T_a=21000.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + N <=> NO + C',
   fr={'Arrhenius', A=2.9e11,  n=0.5, T_a=53630.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO + CO <=> CO2 + N',
   fr={'Arrhenius', A=4.6e8,  n=0.5, T_a=12070.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO + NO <=> N2 + O2',
   fr={'Arrhenius', A=3.07e11,  n=0.0, T_a=33660.0},
   ec={model='from CEA curves',iT=0}
}

-- 8. Associative ionisation

reaction{
   'CO+ + e- <=> C + O',
   fr={'Park', A=1.86e18, n=-0.48, T_a=0.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N + O <=> NO+ + e-',
   fr={'Arrhenius', A=8.8e8, n=1.0, T_a=31900.0},
   ec={model='from CEA curves',iT=0}
   -- New reaction: from Park (1994) scheme
}

reaction{
   'O + O <=> O2+ + e-',
   fr={'Arrhenius', A=7.1e2, n=2.7, T_a=80600.0},
   ec={model='from CEA curves',iT=0}
   -- New reaction: from Park (1994) scheme
}

-- 9. Electron impact ionisation 

reaction{
   'C + e- <=> C+ + e- + e-',
   fr={'Park', A=5.85e27, n=-2.074, T_a=127510.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=127510.0} },
   ec={model='from CEA curves',iT=-1,species='e_minus',mode='translation'}
}

reaction{
   'N + e- <=> N+ + e- + e-',
   fr={'Park', A=1.93e31, n=-2.856, T_a=168970.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=168970.0} },
   ec={model='from CEA curves',iT=-1,species='e_minus',mode='translation'}
}

-- NOTE: the A parameter in the following reaction has been corrected from 8.25e37 to 8.25e27, based on comparison with the park rate

reaction{
   'O + e- <=> O+ + e- + e-',
   fr={'Park', A=8.25e27, n=-2.237, T_a=157840.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=157840.0} },
   ec={model='from CEA curves',iT=-1,species='e_minus',mode='translation'}
}

-- 10. Charge exchange reactions

reaction{
   'CO + C+ <=> CO+ + C',
   fr={'Arrhenius', A=1.0e13, n=0.0, T_a=31400.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO+ + O <=> CO + O+',
   fr={'Arrhenius', A=8.43e13, n=0.0, T_a=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N+ + O <=> N + O+',
   fr={'Arrhenius', A=6.02e11, n=0.0, T_a=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O+ + N <=> N+ + O',
   fr={'Arrhenius', A=7.83e13, n=0.0, T_a=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO+ + C <=> NO + C+',
   fr={'Arrhenius', A=1.0e13,  n=0.0, T_a=23200.0},
   ec={model='from CEA curves',iT=0}
   -- New reaction: from Park (1994) scheme
}

reaction{
   'O2+ + O <=> O+ + O2',
   fr={'Arrhenius', A=4.0e12,  n=-0.09, T_a=18000.0},
   ec={model='from CEA curves',iT=0}
   -- New reaction: from Park (1994) scheme
}

reaction{
   'NO+ + N <=> O+ + N2',
   fr={'Arrhenius', A=3.4e13,  n=-1.08, T_a=12800.0},
   ec={model='from CEA curves',iT=0}
   -- New reaction: from Park (1994) scheme
}

reaction{
   'NO+ + O <=> O2+ + N',
   fr={'Arrhenius', A=7.2e12,  n=0.29, T_a=48600.0},
   ec={model='from CEA curves',iT=0}
   -- New reaction: from Park (1994) scheme
}

reaction{
   'O2 + C+ <=> O2+ + C',
   fr={'Arrhenius', A=1.0e13,  n=0.0, T_a=9400.0},
   ec={model='from CEA curves',iT=0}
   -- New reaction: from Park (1994) scheme
}

-- 11. Ion exchange reactions

reaction{
   'N+ + O2 <=> O+ + NO',
   fr={'Arrhenius', A=1.69e13, n=0.0, T_a=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N+ + NO <=> N2 + O+',
   fr={'Arrhenius', A=6.02e11, n=0.0, T_a=0.0},
   ec={model='from CEA curves',iT=0}
}
