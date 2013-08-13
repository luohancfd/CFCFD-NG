-- Abe1980_26_sp_50r
--
-- 
-- Author: Jeremy Mora
-- Date: 10-May-2012
-- Place: EPFL, Lausanne, Switzerland
--
-- Species:
-- 'C','C2','C3','CO2','CO','CN','CO_plus','C_plus','C2H2','C2H','CH','H','H2','HCN','CHO','N2','N','NO','N_plus','N2_plus','NO_plus','O2','O','O_plus','O2_plus','e_minus' 

scheme_t = {
    update = "chemical kinetic ODE MC",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

-- Dissociation reactions

reaction{
   'C2 + M <=> C + C + M',
   fr={'Park', A=3.7e14, n=0.0, T_a=69900.0, p_name='C2', p_mode='vibration', s_p=0.3, q_name='C2', q_mode='translation'},
efficiencies={C=1.0,C2=1.0,C3=1.0,CO2=0.0,CO=1.0,CN=1.0,CO_plus=0.0,C_plus=1.0,C2H2=0.0,C2H=1.0,CH=0.0,H=1.0,H2=1.0,HCN=0.0,CHO=0.0,N2=1.0,N=1.0,NO=1.0,N_plus=1.0,N2_plus=1.0,NO_plus=1.0,O2=1.0,O=1.0,O_plus=1.0,O2_plus=0.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'C3 + M <=> C + C2 + M',
   fr={'Park', A=1.6e16, n=1.0, T_a=87480.0, p_name='C3', p_mode='vibration', s_p=0.3, q_name='C3', q_mode='translation'},
efficiencies={C=1.0,C2=1.0,C3=1.0,CO2=1.0,CO=1.0,CN=1.0,CO_plus=0.0,C_plus=0.0,C2H2=0.0,C2H=0.0,CH=0.0,H=1.0,H2=0.0,HCN=0.0,CHO=0.0,N2=1.0,N=1.0,NO=1.0,N_plus=0.0,N2_plus=0.0,NO_plus=1.0,O2=0.0,O=1.0,O_plus=0.0,O2_plus=0.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'C2H2 + M <=> C2H + H + M',
   fr={'Park', A=6.96e39, n=-6.1, T_a=67130.0, p_name='C2H2', p_mode='vibration', s_p=0.3, q_name='C2H2', q_mode='translation'},
efficiencies={C=1.0,C2=1.0,C3=1.0,CO2=1.0,CO=1.0,CN=1.0,CO_plus=1.0,C_plus=1.0,C2H2=1.0,C2H=1.0,CH=1.0,H=1.0,H2=1.0,HCN=1.0,CHO=1.0,N2=1.0,N=1.0,NO=1.0,N_plus=1.0,N2_plus=1.0,NO_plus=1.0,O2=1.0,O=1.0,O_plus=1.0,O2_plus=1.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CN + M <=> C + N + M',
   fr={'Park', A=2.5e14, n=0.0, T_a=87740.0, p_name='CN', p_mode='vibration', s_p=0.3, q_name='CN', q_mode='translation'},
efficiencies={C=1.0,C2=1.0,C3=1.0,CO2=0.0,CO=1.0,CN=1.0,CO_plus=0.0,C_plus=1.0,C2H2=0.0,C2H=1.0,CH=0.0,H=1.0,H2=1.0,HCN=0.0,CHO=0.0,N2=1.0,N=1.0,NO=1.0,N_plus=1.0,N2_plus=1.0,NO_plus=1.0,O2=1.0,O=1.0,O_plus=1.0,O2_plus=0.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + M <=> C + O + M',
   fr={'Park', A=2.3e20, n=-1.0, T_a=129000.0, p_name='CO', p_mode='vibration', s_p=0.3, q_name='CO', q_mode='translation'},
efficiencies={C=1.478,C2=1.0,C3=0.0,CO2=1.0,CO=1.0,CN=1.0,CO_plus=0.0,C_plus=0.0,C2H2=0.0,C2H=0.0,CH=0.0,H=0.0,H2=0.0,HCN=0.0,CHO=0.0,N2=1.0,N=1.478,NO=1.0,N_plus=0.0,N2_plus=0.0,NO_plus=0.0,O2=1.0,O=1.478,O_plus=0.0,O2_plus=0.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO2 + M <=> CO + O + M',
   fr={'Park', A=3.5e14, n=0.0, T_a=52525.0, p_name='CO2', p_mode='vibration', s_p=0.3, q_name='CO2', q_mode='translation'},
efficiencies={C=1.0,C2=1.0,C3=1.0,CO2=1.0,CO=1.0,CN=1.0,CO_plus=1.0,C_plus=1.0,C2H2=1.0,C2H=1.0,CH=1.0,H=1.0,H2=1.0,HCN=1.0,CHO=1.0,N2=1.0,N=1.0,NO=1.0,N_plus=1.0,N2_plus=1.0,NO_plus=1.0,O2=1.0,O=1.0,O_plus=1.0,O2_plus=1.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'H2 + M <=> H + H + M',
   fr={'Park', A=2.2e14, n=0.0, T_a=48300.0, p_name='H2', p_mode='vibration', s_p=0.3, q_name='H2', q_mode='translation'},
efficiencies={C=1.0,C2=1.0,C3=1.0,CO2=0.0,CO=1.0,CN=1.0,CO_plus=0.0,C_plus=1.0,C2H2=0.0,C2H=1.0,CH=0.0,H=1.0,H2=2.5,HCN=0.0,CHO=0.0,N2=1.0,N=1.0,NO=1.0,N_plus=1.0,N2_plus=1.0,NO_plus=1.0,O2=1.0,O=1.0,O_plus=1.0,O2_plus=0.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CHO + M <=> H + CO + M',
   fr={'Park', A=1.87e17, n=-1.0, T_a=8560.5, p_name='CHO', p_mode='vibration', s_p=0.3, q_name='CHO', q_mode='translation'},
   efficiencies={C=1.0,C2=1.0,C3=1.0,CO2=1.0,CO=1.0,CN=1.0,CO_plus=1.0,C_plus=1.0,C2H2=1.0,C2H=1.0,CH=1.0,H=1.0,H2=1.0,HCN=1.0,CHO=1.0,N2=1.0,N=1.0,NO=1.0,N_plus=1.0,N2_plus=1.0,NO_plus=1.0,O2=1.0,O=1.0,O_plus=1.0,O2_plus=1.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + M <=> N + N + M',
   fr={'Park', A=7.0e21, n=-1.60, T_a=113200.0, p_name='N2', p_mode='vibration', s_p=0.3, q_name='N2', q_mode='translation'},
   efficiencies={C=4.286,C2=1.0,C3=1.0,CO2=0.0,CO=1.0,CN=1.0,CO_plus=0.0,C_plus=1.0,C2H2=0.0,C2H=1.0,CH=1.0,H=4.286,H2=1.0,HCN=0.0,CHO=0.0,N2=1.0,N=4.286,NO=1.0,N_plus=1.0,N2_plus=1.0,NO_plus=1.0,O2=1.0,O=4.286,O_plus=1.0,O2_plus=0.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO + M <=> N + O + M',
   fr={'Park', A=5.0e15, n=0.0, T_a=75500.0, p_name='NO', p_mode='vibration', s_p=0.3, q_name='NO', q_mode='translation'},
   efficiencies={C=22.0,C2=1.0,C3=0.0,CO2=22.0,CO=1.0,CN=1.0,CO_plus=0.0,C_plus=0.0,C2H2=0.0,C2H=0.0,CH=0.0,H=22.0,H2=0.0,HCN=0.0,CHO=0.0,N2=1.0,N=22.0,NO=22.0,N_plus=0.0,N2_plus=0.0,NO_plus=0.0,O2=1.0,O=22.0,O_plus=0.0,O2_plus=0.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O2 + M <=> O + O + M',
   fr={'Park', A=2.0e21, n=-1.5, T_a=59360.0, p_name='O2', p_mode='vibration', s_p=0.3, q_name='O2', q_mode='translation'},
efficiencies={C=5.0,C2=1.0,C3=1.0,CO2=0.0,CO=1.0,CN=1.0,CO_plus=0.0,C_plus=1.0,C2H2=0.0,C2H=1.0,CH=1.0,H=5.0,H2=1.0,HCN=0.0,CHO=0.0,N2=1.0,N=5.0,NO=1.0,N_plus=1.0,N2_plus=1.0,NO_plus=1.0,O2=1.0,O=5.0,O_plus=1.0,O2_plus=0.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

-- Neutral exchange reactions

reaction{
   'C2 + CO <=> C3 + O',
   fr={'Arrhenius', A=1.2e13, n=0.0, T_a=43240.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'C3 + C <=> C2 + C2',
   fr={'Arrhenius', A=1.0e12, n=0.0, T_a=16400.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CH + O2 <=> CHO + O',
   fr={'Arrhenius', A=6.71e13, n=0.0, T_a=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'C2H + O <=> CH + CO',
   fr={'Arrhenius', A=5.0e13, n=0.0, T_a=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'C2H + O2 <=> CHO + CO',
   fr={'Arrhenius', A=1.0e13, n=0.0, T_a=3524.9},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'C2H2 + H <=> C2H + H2',
   fr={'Arrhenius', A=6.62e13, n=0.0, T_a=14000.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'C2H2 + O2 <=> CHO + CHO',
   fr={'Arrhenius', A=3.98e12, n=0.0, T_a=14106.7},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CN + H2 <=> HCN + H',
   fr={'Arrhenius', A=2.95e5, n=2.5, T_a=1130.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CN + O <=> C + NO',
   fr={'Arrhenius', A=1.6e13, n=0.1, T_a=14600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + C <=> C2 + O',
   fr={'Arrhenius', A=2.0e17, n=-1.0, T_a=58000.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + CO <=> C2 + O2',
   fr={'Arrhenius', A=9.2e11, n=0.75, T_a=163300.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + CO <=> CO2 + C',
   fr={'Arrhenius', A=1.0e3, n=2.0, T_a=72390.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + N <=> C + NO',
   fr={'Arrhenius', A=9.0e16, n=-1.0, T_a=53200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + N <=> CN + O',
   fr={'Arrhenius', A=1.0e14, n=0.0, T_a=38600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + NO <=> CO2 + N',
   fr={'Arrhenius', A=1.0e3, n=2.0, T_a=20980.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + O <=> O2 + C',
   fr={'Arrhenius', A=3.9e13, n=-0.18, T_a=69200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO2 + N <=> CN + O2',
   fr={'Arrhenius', A=3.0e8, n=1.0, T_a=49560.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO2 + O <=> CO + O2',
   fr={'Arrhenius', A=2.1e13, n=0.0, T_a=27800.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CHO + H <=> CO + H2',
   fr={'Arrhenius', A=7.34e13, n=0.0, T_a=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + C <=> CN + N',
   fr={'Arrhenius', A=1.1e14, n=-0.11, T_a=23200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + CO <=> CN + NO',
   fr={'Arrhenius', A=1.0e3, n=2.0, T_a=92010.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + O <=> NO + N',
   fr={'Arrhenius', A=5.7e12, n=0.42, T_a=42938.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO + O <=> O2 + N',
   fr={'Arrhenius', A=8.4e12, n=0.0, T_a=19400.0},
   ec={model='from CEA curves',iT=0}
}

-- Associative ionization reactions

reaction{
   'C + O <=> CO+ + e-',
   fr={'Arrhenius', A=8.8e8, n=1.0, T_a=33100.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N + N <=> N2+ + e-',
   fr={'Arrhenius', A=4.4e7, n=1.5, T_a=67500.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N + O <=> NO+ + e-',
   fr={'Arrhenius', A=5.3e12, n=0.0, T_a=31900.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O + O <=> O2+ + e-',
   fr={'Arrhenius', A=7.1e2, n=2.7, T_a=80600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O2 + N2 <=> NO + NO+ + e-',
   fr={'Arrhenius', A=1.38e20, n=-1.84, T_a=141000.0},
   ec={model='from CEA curves',iT=0}
}

-- Electron-impact ionization reactions

reaction{
   'C + e- <=> C+ + e- + e-',
   fr={'Park', A=3.7e31, n=-3.00, T_a=130720.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=130720.0*2} },
   ec={model='from CEA curves',iT=1}
}

reaction{
   'N + e- <=> N+ + e- + e-',
   fr={'Park', A=2.5e34, n=-3.82, T_a=168200.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=168200.0*2} },
   ec={model='from CEA curves',iT=1}
}

reaction{
   'O + e- <=> O+ + e- + e-',
   fr={'Park', A=3.9e33, n=-3.78, T_a=158500.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=158500.0*2} },
   ec={model='from CEA curves',iT=1}
}

--Charge exchange reactions

reaction{
   'N2 + N+ <=> N + N2+',
   fr={'Arrhenius', A=1.0e12, n=0.5, T_a=12200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + O+ <=> O + N2+',
   fr={'Arrhenius', A=3.4e19, n=-2.0, T_a=23000.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N + NO+ <=> NO + N+',
   fr={'Arrhenius', A=1.0e19, n=-0.93, T_a=61000.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO + M <=> NO+ + e- + M',
   fr={'Arrhenius', A=2.2e15, n=-0.35, T_a=108000.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O + NO+ <=> NO + O+',
   fr={'Arrhenius', A=3.63e15, n=-0.6, T_a=50800.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O + NO+ <=> O2 + N+',
   fr={'Arrhenius', A=1.0e12, n=0.5, T_a=77200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O2 + NO+ <=> NO + O2+',
   fr={'Arrhenius', A=2.4e13, n=0.41, T_a=32600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O + O2+ <=> O+ + O2',
   fr={'Arrhenius', A=4.0e12, n=-0.09, T_a=18000.0},
   ec={model='from CEA curves',iT=0}
}
