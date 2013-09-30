-- Park2001-NEQ-20sp.lua
--
-- Park, C.S. (2001)
-- Chemical-Kinetic Parameters of Hyperbolic Earth Entry
-- JTHT Volume 15 Number 1 pp 76-90 Jan-Mar. 2001
--
-- Author: Ojas Joshi
-- Date: 06-Sep-2011
-- Place: EPFL, Lausanne, Switzerland
--
-- Species:
-- 'C', 'O', 'N', 'H', 'CO', 'C2', 'N2', 'CN', 'NO', 'O2', 'H2', 'C3', 'C2H', 'C_plus', 'O_plus', 'H_plus', 'N_plus', 'NO_plus', 'N2_plus', 'e_minus' 

scheme_t = {
    update = "chemical kinetic ODE MC",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

reaction{
   'O2 + M <=> O + O + M',
   fr={'Park', A=2.0e21, n=-1.5, T_a=59360.0, p_name='O2', p_mode='vibration', s_p=0.3, q_name='O2', q_mode='translation'},
efficiencies={C=5.0,O=5.0,CO=1.0,C2=1.0,O2=1.0,C3=1.0,C_plus=1.0,O_plus=1.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'C2 + M <=> C + C + M',
   fr={'Park', A=3.7e14, n=0.0, T_a=69900.0, p_name='C2', p_mode='vibration', s_p=0.3, q_name='C2', q_mode='translation'},
efficiencies={C=1.0,O=1.0,CO=1.0,C2=1.0,O2=1.0,C3=1.0,C_plus=1.0,O_plus=1.0,e_minus=0.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + O <=> O2 + C',
   fr={'Arrhenius', A=3.9e13, n=-0.18, T_a=69200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + C <=> C2 + O',
   fr={'Arrhenius', A=2.0e17, n=-1, T_a=58000.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'CO + C2 <=> C3 + O',
   fr={'Arrhenius', A=1.0e12, n=0.0, T_a=41200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'C3 + C <=> C2 + C2',
   fr={'Arrhenius', A=1.0e12, n=0.0, T_a=16400.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O + e- <=> O+ + e- + e-',
   fr={'Park', A=3.9e33, n=-3.78, T_a=158500.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=158500.0*2} },
   ec={model='from CEA curves',iT=1}
}

reaction{
   'C + e- <=> C+ + e- + e-',
   fr={'Park', A=3.7e31, n=-3.00, T_a=130720.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=130720.0*2} },
   ec={model='from CEA curves',iT=1}
}

