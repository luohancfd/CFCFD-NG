-- Park93-1T-5sp.lua
--
-- Single temperature, neutral reactions from:
--
-- Park, C.S. (1993)
-- Review of Chemical-Kinetic problems of future NASA missions, II: Earth Entries
-- JTHT Volume 7 Number 3 pp 385-398 July-Sept. 1993
--
-- Author: Daniel F. Potter
-- Date: 01-Dec-2009
-- Place: UQ, Brisbane, Australia
--
-- History: 
--   01-Dec-2009: - Port from gas_models2 input file, correcting effeciencies
--   02-Dec-2009: - Ion reactions
--   07-Dec-2009: - Nonequilibrium rates
--   07-Aug-2013: - Knab CVCV model for dissociation and neutral exchange reactions
--   11-Spe-2013: - Using alpha and U values from Kanne et al

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
   'N2 + M <=> N + N + M',
   fr={'Knab_et_al', A=7.0e21, n=-1.60, T_a=113200.0, v_name='N2', U=113200.0/6., alpha=0.9},
   chemistry_energy_coupling={ {species='N2', mode='vibration', model='Knab_et_al', A=113200, U=113200.0/6., alpha=0.9} },
   efficiencies={N2=1.0,N2_plus=1.0,O2=1.0,O2_plus=1.0,NO=1.0,NO_plus=1.0,N=4.2857,N_plus=4.2857,O=4.2857,O_plus=4.2857,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + e- <=> N + N + e-',
   fr={'Park', A=1.0e25, n=-1.60, T_a=113200.0, p_name='N2', p_mode='vibration', s_p=0.5, q_name='e_minus', q_mode='translation'},
   ec={model='from CEA curves',iT=1}
}

reaction{
   'O2 + M <=> O + O + M',
   fr={'Knab_et_al', A=2.0e21, n=-1.5, T_a=59500.0, v_name='O2', U=59500.0/6., alpha=0.8},
   chemistry_energy_coupling={ {species='O2', mode='vibration', model='Knab_et_al', A=59500.0, U=59500.0/6., alpha=0.8} },
   efficiencies={N2=1.0,N2_plus=1.0,O2=1.0,O2_plus=1.0,NO=1.0,NO_plus=1.0,N=5.0,N_plus=5.0,O=5.0,O_plus=5.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO + M <=> N + O + M',
   fr={'Knab_et_al', A=5.0e15, n=0.0, T_a=75500.0, v_name='NO', U=75500.0/5., alpha=0.8},
   chemistry_energy_coupling={ {species='NO', mode='vibration', model='Knab_et_al', A=75500.0, U=75500.0/5., alpha=0.8} },
   efficiencies={N2=1.0,N2_plus=1.0,O2=1.0,O2_plus=1.0,NO=1.0,NO_plus=1.0,N=20.0,N_plus=20.0,O=20.0,O_plus=20.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

-- Neutral exchange reactions

reaction{
   'NO + O <=> O2 + N',
   fr={'Knab_et_al', A=8.4e12, n=0.0, T_a=19450.0, v_name='NO', U=75500.0/5., alpha=0.1},
   chemistry_energy_coupling={ {species='NO', mode='vibration', model='Knab_et_al', A=19450.0, U=75500.0/5., alpha=0.1},
                               {species='O2', mode='vibration', model='Knab_et_al', A=19450.0, U=59500.0/5., alpha=0.1} },
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + O <=> NO + N',
   fr={'Knab_et_al', A=6.4e17, n=-1.0, T_a=38400.0, v_name='N2', U=113200.0/5., alpha=0.1},
   chemistry_energy_coupling={ {species='NO', mode='vibration', model='Knab_et_al', A=38400.0,  U=75500.0/5., alpha=0.1},
                               {species='N2', mode='vibration', model='Knab_et_al', A=38400.0, U=113200.0/5., alpha=0.1} },
   ec={model='from CEA curves',iT=0}
}

-- Associative ionization reactions

reaction{
   'N + O <=> NO+ + e-',
   fr={'Arrhenius', A=8.8e8, n=1.0, T_a=31900.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O + O <=> O2+ + e-',
   fr={'Arrhenius', A=7.1e2, n=2.7, T_a=80600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N + N <=> N2+ + e-',
   fr={'Arrhenius', A=4.4e7, n=1.5, T_a=67500.0},
   ec={model='from CEA curves',iT=0}
}

-- Charge exchange reactions

reaction{
   'NO+ + O <=> N+ + O2',
   fr={'Arrhenius', A=1.0e12, n=0.5, T_a=77200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N+ + N2 <=> N2+ + N',
   fr={'Arrhenius', A=1.0e12, n=0.5, T_a=12200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O2+ + N <=> N+ + O2',
   fr={'Arrhenius', A=8.7e13, n=0.14, T_a=28600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O+ + NO <=> N+ + O2',
   fr={'Arrhenius', A=1.4e5, n=1.9, T_a=26600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O2+ + N2 <=> N2+ + O2',
   fr={'Arrhenius', A=9.9e12, n=0.0, T_a=40700.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O2+ + O <=> O+ + O2',
   fr={'Arrhenius', A=4.0e12, n=-0.09, T_a=18000.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO+ + N <=> O+ + N2',
   fr={'Arrhenius', A=3.4e13, n=-1.08, T_a=12800.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO+ + O2 <=> O2+ + NO',
   fr={'Arrhenius', A=2.4e13,  n=0.41, T_a=32600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO+ + O <=> O2+ + N',
   fr={'Arrhenius', A=7.2e12,  n=0.29, T_a=48600.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O+ + N2 <=> N2+ + O',
   fr={'Arrhenius', A=9.1e11,  n=0.36, T_a=22800.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO+ + N <=> N2+ + O',
   fr={'Arrhenius', A=7.2e13,  n=0.00, T_a=35500.0},
   ec={model='from CEA curves',iT=0}
}

-- Electron-impact ionization reactions

reaction{
   'O + e- <=> O+ + e- + e-',
   fr={'Park', A=3.9e33, n=-3.78, T_a=158500.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=51720.0} },
   ec={model='from CEA curves',iT=1}
}

reaction{
   'N + e- <=> N+ + e- + e-',
   fr={'Park', A=2.5e34, n=-3.82, T_a=168600.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=48710.0} },
   ec={model='from CEA curves',iT=1}
}

-- Radiative recombination reactions [omitted]

-- reaction{
--     'O+ + e- <=> O',
--     fr={'Park', A=1.07e11, n=-0.52, T_a=0.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
--     ec={model='from CEA curves',iT=1}
-- }

-- reaction{
--     'N+ + e- <=> N',
--     fr={'Park', A=1.52e11, n=-0.48, T_a=0.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
--     ec={model='from CEA curves',iT=1}
-- }
