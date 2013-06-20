-- Park2001-20sp-24r.lua
--
-- Park, C.S. (2001)
-- Chemical-Kinetic Parameters of Hyperbolic Earth Entry
-- JTHT Volume 15 Number 1 pp 76-90 Jan-Mar. 2001
--
-- Author: Ojas Joshi
-- Date: 06-Sep-2011
-- Place: EPFL, Lausanne, Switzerland
-- Updated by Elise Fahy
-- Feb 2013, EPFL, Lausanne, Switzerland
--
-- Species:
-- 'C', 'O', 'N', 'H', 'CO', 'C2', 'N2', 'CN', 'NO', 'O2', 'H2', 'C3', 'C2H', 'C_plus', '
-- O_plus', 'H_plus', 'N_plus', 'NO_plus', 'N2_plus', 'e_minus' 

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
   fr={'Park', A=7.0e21, n=-1.60, T_a=113200.0, p_name='N2', p_mode='vibration', s_p=0.3, q_name='N2', q_mode='translation'},
   efficiencies={C=4.2857,O=4.2857,N=4.2857,H=4.2857,CO=1.0,C2=1.0,N2=1.0,CN=1.0,NO=1.0,O2=1.0,H2=1.0,C3=1.0,C2H=1.0,C_plus=1.0,
   H_plus=1.0,O_plus=1.0,N_plus=1.0,NO_plus=1.0,N2_plus=1.0,e_minus=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2 + e- <=> N + N + e-',
   fr={'Park', A=3.0e24, n=-1.60, T_a=113200.0, p_name='N2', p_mode='vibration', s_p=0.3, q_name='e_minus', q_mode='translation'},
   ec={model='from thermo',iT=0}
}

reaction{
   'O2 + M <=> O + O + M',
   fr={'Park', A=2.0e21, n=-1.5, T_a=59360.0, p_name='O2', p_mode='vibration', s_p=0.3, q_name='O2', q_mode='translation'},
   efficiencies={C=5.0, O=5.0, N=5.0, H=5.0, CO=1.0, C2=1.0, N2=1.0, CN=1.0, NO=1.0, O2=1.0, H2=1.0, C3=1.0, C2H=1.0, C_plus=1.0,
   H_plus=1.0, O_plus=1.0, N_plus=1.0, NO_plus=1.0, N2_plus=1.0, e_minus=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'C2 + M <=> C + C + M',
   fr={'Park', A=3.7e14, n=0.0, T_a=69900.0, p_name='C2', p_mode='vibration', s_p=0.3, q_name='C2', q_mode='translation'},
   efficiencies={C=1.0, O=1.0, N=1.0, H=1.0, CO=1.0, C2=1.0, N2=1.0, CN=1.0, NO=1.0, O2=1.0, H2=1.0, C3=1.0, C2H=1.0, C_plus=1.0, 
   H_plus=1.0, O_plus=1.0, N_plus=1.0, NO_plus=1.0, N2_plus=1.0, e_minus=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CN + M <=> C + N + M',
   fr={'Park', A=2.5e14, n=0.0, T_a=87740.0, p_name='CN', p_mode='vibration', s_p=0.3, q_name='CN', q_mode='translation'},
   efficiencies= {C=1.0, O=1.0, N=1.0, H=1.0, CO=1.0, C2=1.0, N2=1.0, CN=1.0, NO=1.0, O2=1.0, H2=1.0, C3=1.0, C2H=1.0, C_plus=1.0, 
   H_plus=1.0, O_plus=1.0, N_plus=1.0, NO_plus=1.0, N2_plus=1.0, e_minus=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'H2 + M <=> H + H + M',
   fr={'Park', A=2.2e14, n=0.0, T_a=48300.0, p_name='H2', p_mode='vibration', s_p=0.3, q_name='H2', q_mode='translation'},
   efficiencies={C=1.0, O=1.0, N=1.0, H=1.0, CO=1.0, C2=1.0, N2=1.0, CN=1.0, NO=1.0, O2=1.0, H2=1.0, C3=1.0, C2H=1.0, C_plus=1.0, 
   H_plus=1.0, O_plus=1.0, N_plus=1.0, NO_plus=1.0, N2_plus=1.0, e_minus=0.0},
   ec={model='from thermo',iT=0}
}

-- Neutral exchange reactions

reaction{
   'N2 + O <=> NO + N',
   fr={'Arrhenius', A=5.7e12, n=0.42, T_a=42938.0},
  ec={model='from thermo',iT=0}
}

reaction{
   'NO + O <=> O2 + N',
   fr={'Arrhenius', A=8.4e12, n=0.0, T_a=19400.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO + C <=> C2 + O',
   fr={'Arrhenius', A=2.0e17, n=-1, T_a=58000.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO + O <=> O2 + C',
   fr={'Arrhenius', A=3.9e13, n=-0.18, T_a=69200.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO + N <=> CN + O',
   fr={'Arrhenius', A=1.0e14, n=0.0, T_a=38600.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2 + C <=> CN + N',
   fr={'Arrhenius', A=1.1e14, n=-0.11, T_a=23200.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CN + O <=> NO + C',
   fr={'Arrhenius', A=1.6e13, n=0.1, T_a=14600.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CN + C <=> C2 + N',
   fr={'Arrhenius', A=5.0e13, n=0.0, T_a=13000.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO + C2 <=> C3 + O',
   fr={'Arrhenius', A=1.0e12, n=0.0, T_a=41200.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'C3 + N <=> CN + C2',
   fr={'Arrhenius', A=1.0e12, n=0.0, T_a=34200.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'C3 + C <=> C2 + C2',
   fr={'Arrhenius', A=1.0e12, n=0.0, T_a=16400.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'C2H + H <=> C2 + H2',
   fr={'Arrhenius', A=1.0e12, n=0.0, T_a=16770.0},
   ec={model='from thermo',iT=0}
}

-- Associative ionization reactions

reaction{
   'N + O <=> NO+ + e-',
   fr={'Arrhenius', A=5.3e12, n=0.0, T_a=31900.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N + N <=> N2+ + e-',
   fr={'Arrhenius', A=4.4e7, n=1.5, T_a=67500.0},
   ec={model='from thermo',iT=0}
}

-- Electron-impact ionization reactions

reaction{
   'O + e- <=> O+ + e- + e-',
   fr={'Park', A=3.9e33, n=-3.78, T_a=158500.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=158500.0*2} },
   ec={model='from thermo',iT=0}
}

reaction{
   'N + e- <=> N+ + e- + e-',
   fr={'Park', A=2.5e34, n=-3.82, T_a=168200.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=168200.0*2} },
   ec={model='from thermo',iT=0}
}

reaction{
   'C + e- <=> C+ + e- + e-',
   fr={'Park', A=3.7e31, n=-3.00, T_a=130720.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=130720.0*2} },
   ec={model='from thermo',iT=0}
}

reaction{
   'H + e- <=> H+ + e- + e-',
   fr={'Park', A=2.2e30, n=-2.80, T_a=157800.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=157800.0*2} },
   ec={model='from thermo',iT=0}
}
