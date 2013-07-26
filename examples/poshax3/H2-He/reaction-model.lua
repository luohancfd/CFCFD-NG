scheme_t = {
    update = "chemical kinetic ODE MC",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

Na = 6.022e23

-- H2 dissociation
-- Ref: JTHT Vol 26 No 2 April-June 2012 (Park: Nonequilibrium ionization and radiation in Hydrogen-Helium mixtures)

reaction{
   'H2 + H2 <=> H + H + H2',
   fr={'Park', A=2.967e-7*Na, n=-0.5165, T_a=52530.0, p_name='H2', p_mode='translation', s_p=0.5, q_name='H2', q_mode='translation' },
   ec={model='from CEA curves',iT=0}
}

reaction{
   'H2 + H <=> H + H + H',
   fr={'Park', A=3.18e-4*Na, n=-1.0735, T_a=55105.0, p_name='H2', p_mode='translation', s_p=0.5, q_name='H', q_mode='translation' },
   ec={model='from CEA curves',iT=0}
}

reaction{
   'H2 + He <=> H + H + He',
   -- for sqrt(T.Tv) less than 10000 K
   -- fr={'Park', A=9.2719e-4*Na, n=-1.4773, T_a=55105.0, p_name='H2', p_mode='translation', s_p=0.5, q_name='He', q_mode='translation' },
   -- for sqrt(T.Tv) greater than 10000 K
   fr={'Park', A=1.4568e-7*Na, n=-0.5018, T_a=55105.0, p_name='H2', p_mode='translation', s_p=0.5, q_name='He', q_mode='translation' },
   ec={model='from CEA curves',iT=0}
}

-- H ionisation
-- Ref: JTHT Vol 15 No 1 2001 (Park et al: Chemical-Kinetic parameters of hyperbolic Earth entry)

reaction{
   'H + e- <=> H+ + e- + e-',
   fr={'Park', A=2.20e30, n=-2.8, T_a=157800.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=157800.0} },
   ec={model='from CEA curves',iT=-1,species='e_minus',mode='translation'}
}
