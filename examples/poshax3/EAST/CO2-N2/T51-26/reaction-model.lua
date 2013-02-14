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
--   06-October-2012: added reactions from Park involving NO+ and O2+
--   07-October-2012: adding radiative transitions
--   07-October-2012: taking out associative ionization electron coupling

scheme_t = {
    update = "chemical kinetic ODE MC",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

-- 1. Heavy particle-impact electronic excitation

-- 1.1 CN

-- 1.1.1 CN(X) + M -> CN(A) + M

reaction{
    'CN_X + N_4So <=> CN_A + N_4So',
    fr={'Park', A=3.30e15, n=0.047, T_a=18988.0, p_name='CN_X', p_mode='vibration', s_p=0.5, q_name='CN_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_X + O_3P <=> CN_A + O_3P',
    fr={'Park', A=2.71e15, n=0.062, T_a=18794.0, p_name='CN_X', p_mode='vibration', s_p=0.5, q_name='CN_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_X + N2_X <=> CN_A + N2_X',
    fr={'Park', A=1.76e24, n=-1.899, T_a=42648.0, p_name='CN_X', p_mode='vibration', s_p=0.5, q_name='CN_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_X + O2 <=> CN_A + O2',
    fr={'Park', A=1.08e15, n=0.132, T_a=17943.0, p_name='CN_X', p_mode='vibration', s_p=0.5, q_name='CN_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_X + NO <=> CN_A + NO',
    fr={'Park', A=1.16e15, n=0.127, T_a=18010.0, p_name='CN_X', p_mode='vibration', s_p=0.5, q_name='CN_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_X + CO_X <=> CN_A + CO_X',
    fr={'Park', A=1.44e18, n=-0.555, T_a=26310.0, p_name='CN_X', p_mode='vibration', s_p=0.5, q_name='CN_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

-- 1.1.2 CN(A) + M -> CN(B) + M

reaction{
    'CN_A + N_4So <=> CN_B + N_4So',
    fr={'Park', A=5.10e14, n=0.041, T_a=29098.0, p_name='CN_A', p_mode='vibration', s_p=0.5, q_name='CN_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_A + O_3P <=> CN_B + O_3P',
    fr={'Park', A=4.11e14, n=0.057, T_a=28904, p_name='CN_A', p_mode='vibration', s_p=0.5, q_name='CN_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_A + N2_X <=> CN_B + N2_X',
    fr={'Park', A=3.79e13, n=0.271, T_a=26301.0, p_name='CN_A', p_mode='vibration', s_p=0.5, q_name='CN_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_A + O2 <=> CN_B + O2',
    fr={'Park', A=3.25e13, n=0.283, T_a=26159.0, p_name='CN_A', p_mode='vibration', s_p=0.5, q_name='CN_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_A + NO <=> CN_B + NO',
    fr={'Park', A=1.77e14, n=0.122, T_a=28121.0, p_name='CN_A', p_mode='vibration', s_p=0.5, q_name='CN_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_A + CO_X <=> CN_B + CO_X',
    fr={'Park', A=2.22e17, n=-0.561, T_a=36420.0, p_name='CN_A', p_mode='vibration', s_p=0.5, q_name='CN_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

-- 1.1.2 CN(X) + M -> CN(B) + M

reaction{
    'CN_X + M <=> CN_B + M',
    fr={'Park', A=1.80e11, n=0.5, T_a=37000.0, p_name='CN_X', p_mode='vibration', s_p=0.5, q_name='CN_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

-- 1.2 CO

-- 1.2.1 CO(X) + M -> CO(a) + M

reaction{
    'CO_X + N_4So <=> CO_a3 + N_4So',
    fr={'Park', A=5.40e11, n=0.504, T_a=70287.0, p_name='CO_X', p_mode='vibration', s_p=0.5, q_name='CO_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_X + O_3P <=> CO_a3 + O_3P',
    fr={'Park', A=5.40e11, n=0.504, T_a=70287.0, p_name='CO_X', p_mode='vibration', s_p=0.5, q_name='CO_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_X + N2_X <=> CO_a3 + N2_X',
    fr={'Park', A=9.01e11, n=0.504, T_a=70287.0, p_name='CO_X', p_mode='vibration', s_p=0.5, q_name='CO_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_X + O2 <=> CO_a3 + O2',
    fr={'Park', A=3.41e13, n=0.504, T_a=70287.0, p_name='CO_X', p_mode='vibration', s_p=0.5, q_name='CO_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_X + NO <=> CO_a3 + NO',
    fr={'Park', A=6.14e13, n=0.504, T_a=70287.0, p_name='CO_X', p_mode='vibration', s_p=0.5, q_name='CO_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_X + CO_X <=> CO_a3 + CO_X',
    fr={'Park', A=1.02e13, n=0.504, T_a=70287.0, p_name='CO_X', p_mode='vibration', s_p=0.5, q_name='CO_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

-- 1.2.2 CO(X) + M -> CO(A) + M

-- mod: schofield rate

reaction{
    'CO_X + M <=> CO_A + M',
    fr={'Park',  A=7.61e+13, n=0.344, T_a=93669.0, p_name='CO_X', p_mode='vibration', s_p=0.5, q_name='CO_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

-- 1.2.3 CO(A) + M -> CO(a) + M

-- reaction{
--     'CO_A + M <=> CO_a3 + M',
--     fr={'Arrhenius', A=6.00e13, n=0.0, T_a=0.0},
--     ec={model='from thermo',iT=0}
-- }

-- NOTE: reforming this reaction in the other direction so we can use sqrt(T.Tv) as the rate controlling temperature

reaction{
    'CO_a3 + M <=> CO_A + M',
    fr={'Park',  A=3.00e13, n=0.0, T_a=23580.0, p_name='CO_a3', p_mode='vibration', s_p=0.5, q_name='CO_a3', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

-- 1.3 N2

-- 1.3.1 N2(X) + M -> N2(A) + M

reaction{
    'N2_X + N_4So <=> N2_A + N_4So',
    fr={'Park',  A=3.17e12, n=0.507, T_a=72737.0, p_name='N2_X', p_mode='vibration', s_p=0.5, q_name='N2_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_X + O_3P <=> N2_A + O_3P',
    fr={'Park',  A=4.39e12, n=0.507, T_a=72737.0, p_name='N2_X', p_mode='vibration', s_p=0.5, q_name='N2_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_X + N2_X <=> N2_A + N2_X',
    fr={'Park',  A=9.43e6, n=0.507, T_a=72737.0, p_name='N2_X', p_mode='vibration', s_p=0.5, q_name='N2_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_X + O2 <=> N2_A + O2',
    fr={'Park',  A=6.00e11, n=0.507, T_a=72737.0, p_name='N2_X', p_mode='vibration', s_p=0.5, q_name='N2_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_X + NO <=> N2_A + NO',
    fr={'Park',  A=8.48e12, n=0.507, T_a=72737.0, p_name='N2_X', p_mode='vibration', s_p=0.5, q_name='N2_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_X + CO_X <=> N2_A + CO_X',
    fr={'Park',  A=3.54e11, n=0.507, T_a=72737.0, p_name='N2_X', p_mode='vibration', s_p=0.5, q_name='N2_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

-- 1.3.2 N2(A) + M -> N2(B) + M

reaction{
    'N2_A + N_4So <=> N2_B + N_4So',
    fr={'Park',  A=2.67e15, n=0.016, T_a=19253.0, p_name='N2_A', p_mode='vibration', s_p=0.5, q_name='N2_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_A + O_3P <=> N2_B + O_3P',
    fr={'Park', A=2.16e15, n=0.032, T_a=19054.0, p_name='N2_A', p_mode='vibration', s_p=0.5, q_name='N2_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_A + N2_X <=> N2_B + N2_X',
    fr={'Park', A=8.51e18, n=-0.777, T_a=28892.0, p_name='N2_A', p_mode='vibration', s_p=0.5, q_name='N2_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_A + O2 <=> N2_B + O2',
    fr={'Park', A=1.31e18, n=-0.601, T_a=26760.0, p_name='N2_A', p_mode='vibration', s_p=0.5, q_name='N2_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_A + NO <=> N2_B + NO',
    fr={'Park', A=4.43e13, n=0.388, T_a=14722.0, p_name='N2_A', p_mode='vibration', s_p=0.5, q_name='N2_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_A + CO_X <=> N2_B + CO_X',
    fr={'Park', A=7.46e21, n=-1.427, T_a=36800.0, p_name='N2_A', p_mode='vibration', s_p=0.5, q_name='N2_A', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

-- 1.4 C2

-- 1.4.1 C2(X) + M -> C2(d) + M

reaction{
    'C2_X + M <=> C2_d3 + M',
    fr={'Park', A=5.20e13, n=0.15, T_a=28807.0, p_name='C2_X', p_mode='vibration', s_p=0.5, q_name='C2_X', q_mode='translation'},
    ec={model='from thermo',iT=0}
}

-- 1.5 O

reaction{
    'O_3P + M <=> O_1D + M',
    fr={'Arrhenius', A=3.98e13, n=0.1454, T_a=22722.0},
    ec={model='from thermo',iT=0}
}

reaction{
    'O_3P + M <=> O_1S + M',
    fr={'Arrhenius', A=5.20e13, n=0.1454, T_a=48515.0},
    ec={model='from thermo',iT=0}
}

reaction{
    'O_1D + M <=> O_1S + M',
    fr={'Arrhenius', A=4.16e13, n=0.1454, T_a=25794.0},
    ec={model='from thermo',iT=0}
}

-- 1.6 N

reaction{
    'N_4So + M <=> N_2Do + M',
    fr={'Arrhenius', A=4.44e13, n=0.1454, T_a=27669.0},
    ec={model='from thermo',iT=0}
}

reaction{
    'N_4So + M <=> N_2Po + M',
    fr={'Arrhenius', A=5.12e13, n=0.1454, T_a=41500.0},
    ec={model='from thermo',iT=0}
}

reaction{
    'N_2Do + M <=> N_2Po + M',
    fr={'Arrhenius', A=3.47e13, n=0.1454, T_a=13831.0},
    ec={model='from thermo',iT=0}
}

-- 1.7 C

reaction{
    'C_3P + M <=> C_1D + M',
    fr={'Arrhenius', A=3.72e13, n=0.1454, T_a=14625.0},
    ec={model='from thermo',iT=0}
}

reaction{
    'C_3P + M <=> C_1S + M',
    fr={'Arrhenius', A=4.86e13, n=0.1454, T_a=31109.0},
    ec={model='from thermo',iT=0}
}

reaction{
    'C_3P + M <=> C_5So + M',
    fr={'Arrhenius', A=5.69e13, n=0.1454, T_a=48502.0},
    ec={model='from thermo',iT=0}
}

reaction{
    'C_1D + M <=> C_1S + M',
    fr={'Arrhenius', A=3.88e13, n=0.1454, T_a=16484.0},
    ec={model='from thermo',iT=0}
}

reaction{
    'C_1D + M <=> C_5So + M',
    fr={'Arrhenius', A=5.01e13, n=0.1454, T_a=33878.0},
    ec={model='from thermo',iT=0}
}

reaction{
    'C_1S + M <=> C_5So + M',
    fr={'Arrhenius', A=3.95e13, n=0.1454, T_a=17393.0},
    ec={model='from thermo',iT=0}
}

-- 2. Electron-impact electronic excitation

-- 2.1 CN

reaction{
    'CN_X + e- <=> CN_A + e-',
    fr={'Park', A=1.89e11, n=0.733, T_a=12979.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='CN_X', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'CN_X + e- <=> CN_B + e-',
    fr={'Park', A=3.46e12, n=0.483, T_a=36964.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='CN_X', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'CN_A + e- <=> CN_B + e-',
    fr={'Park', A=1.42e18, n=-0.402, T_a=25345.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='CN_X', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

-- 2.2 CO

reaction{
    'CO_X + e- <=> CO_a3 + e-',
    fr={'Park', A=2.59e14, n=0.242, T_a=70255.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='CO_X', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'CO_X + e- <=> CO_A + e-',
    fr={'Park', A=1.14e14, n=0.341, T_a=93739.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='CO_X', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'CO_a3 + e- <=> CO_A + e-',
    fr={'Park', A=1.55e17, n=-0.401, T_a=25172.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='CO_a3', q_mode='translation'},
    br={'Park', A=4.00e15, n=0.0,    T_a=0.0,     p_name='e_minus', p_mode='translation', s_p=1.0, q_name='CO_A', q_mode='translation'},
}

-- Omitting CO(A) + e- -> CO(a) + e- as this is the reverse of the above reaction

-- 2.3 N2

reaction{
    'N2_X + e- <=> N2_A + e-',
    fr={'Park', A=3.30e13, n=0.55, T_a=57700.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='N2_X', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'N2_X + e- <=> N2_B + e-',
    fr={'Park', A=1.12e14, n=0.50, T_a=75288.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='N2_X', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'N2_X + e- <=> N2_C + e-',
    fr={'Park', A=5.33e20, n=-1.13, T_a=124966.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='N2_X', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

-- 2.3 C2

reaction{
    'C2_X + e- <=> C2_d3 + e-',
    fr={'Park', A=7.82e15, n=0.15, T_a=28807.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='C2_X', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

-- 2.4 O

reaction{
    'O_3P + e- <=> O_1D + e-',
    fr={'Park', A=6.68e14, n=-0.046, T_a=24061.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='O_3P', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'O_3P + e- <=> O_1S + e-',
    fr={'Park', A=9.15e12, n=-0.501, T_a=51148.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='O_3P', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'O_1D + e- <=> O_1S + e-',
    fr={'Park', A=2.05e19, n=-1.06, T_a=29888.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='O_1D', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

-- 2.5 N

reaction{
    'N_4So + e- <=> N_2Do + e-',
    fr={'Park', A=5.04e14, n=-0.00911, T_a=28851.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='N_4So', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'N_4So + e- <=> N_2Po + e-',
    fr={'Park', A=2.70e14, n=0.0708, T_a=42571.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='N_4So', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'N_2Do + e- <=> N_2Po + e-',
    fr={'Park', A=1.13e19, n=-1.085, T_a=17714.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='N_2Do', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

-- 2.6 C

reaction{
    'C_3P + e- <=> C_1D + e-',
    fr={'Park', A=1.17e19, n=-1.085, T_a=17714.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='C_3P', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'C_3P + e- <=> C_1S + e-',
    fr={'Park', A=6.56e17, n=-0.69, T_a=33887.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='C_3P', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'C_3P + e- <=> C_5So + e-',
    fr={'Park', A=1.55e17, n=-0.502, T_a=51067.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='C_3P', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'C_1D + e- <=> C_1S + e-',
    fr={'Park', A=1.26e19, n=-1.075, T_a=19724.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='C_1D', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'C_1D + e- <=> C_5So + e-',
    fr={'Park', A=3.70e14, n=0.0308, T_a=35012.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='C_1D', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
    'C_1S + e- <=> C_5So + e-',
    fr={'Park', A=1.32e19, n=-1.072, T_a=20725.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='C_1S', q_mode='translation'},
    ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

-- 3. Heavy particle-impact dissociation

-- 3.1 C2

reaction{
    'C2_X + M <=> C_3P + C_3P + M',
    fr={'Park', A=1.50e16, n=0.0, T_a=71600.0, p_name='C2_X', p_mode='translation', s_p=0.5, q_name='C2_X', q_mode='vibration'},
    chemistry_energy_coupling={ {species='C2_X', mode='vibration', model='TreanorMarrone', U=71600.0/3.0 } },
    ec={model='from thermo',iT=0}
}

-- 3.2 CO2

reaction{
    'CO2 + M <=> CO_X + O_3P + M',
    fr={'Park', A=1.40e22, n=-1.50, T_a=63275.0, p_name='CO2', p_mode='translation', s_p=0.5, q_name='CO2', q_mode='vibration'},
    efficiencies={CO2=0.0,CO_X=0.0,CO_a3=0.0,CO_A=0.0,CO_plus=0.0,O2=0.0,N2_X=0.0,N2_A=0.0,N2_B=0.0,N2_C=0.0,NO=0.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=1.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=1.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=1.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO2', mode='vibration', model='TreanorMarrone', U=63275.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'CO2 + M <=> CO_X + O_3P + M',
    fr={'Park', A=6.90e21, n=-1.50, T_a=63275.0, p_name='CO2', p_mode='translation', s_p=0.5, q_name='CO2', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=1.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=1.0,CN_A=1.0,CN_B=1.0,C2_X=1.0,C2_d3=1.0,
                  C_3P=0.0,C_1D=0.0,C_1S=0.0,C_5So=0.0,C_plus=0.0,O_3P=0.0,O_1D=0.0,O_1S=0.0,O_plus=0.0,N_4So=0.0,N_2Do=0.0,N_2Po=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO2', mode='vibration', model='TreanorMarrone', U=63275.0/3.0 } },
    ec={model='from thermo',iT=0}
}

-- 3.3 CN

reaction{
    'CN_X + M <=> C_3P + N_4So + M',
    fr={'Park', A=2.53e14, n=0.0, T_a=88672.0, p_name='CN_X', p_mode='translation', s_p=0.5, q_name='CN_X', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=1.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=1.0,CN_A=1.0,CN_B=1.0,C2_X=1.0,C2_d3=1.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=1.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=1.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=1.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CN_X', mode='vibration', model='TreanorMarrone', U=88672.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_A + M <=> C_3P + N_4So + M',
    fr={'Park', A=2.53e14, n=0.0, T_a=75371.0, p_name='CN_A', p_mode='translation', s_p=0.5, q_name='CN_A', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=1.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=1.0,CN_A=1.0,CN_B=1.0,C2_X=1.0,C2_d3=1.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=1.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=1.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=1.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CN_A', mode='vibration', model='TreanorMarrone', U=75371.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_B + M <=> C_3P + N_4So + M',
    fr={'Park', A=2.53e14, n=0.0, T_a=51620.0, p_name='CN_B', p_mode='translation', s_p=0.5, q_name='CN_B', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=1.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=1.0,CN_A=1.0,CN_B=1.0,C2_X=1.0,C2_d3=1.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=1.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=1.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=1.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CN_B', mode='vibration', model='TreanorMarrone', U=51620.0/3.0 } },
    ec={model='from thermo',iT=0}
}

-- 3.3 CO

reaction{
    'CO_X + M <=> C_3P + O_3P + M',
    fr={'Park', A=3.40e20, n=-1.00, T_a=128755.0, p_name='CO_X', p_mode='translation', s_p=0.5, q_name='CO_X', q_mode='vibration'},
    efficiencies={CO2=0.0,CO_X=0.0,CO_a3=0.0,CO_A=0.0,CO_plus=0.0,O2=0.0,N2_X=0.0,N2_A=0.0,N2_B=0.0,N2_C=0.0,NO=0.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=0.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=0.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO_X', mode='vibration', model='TreanorMarrone', U=128755.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_X + M <=> C_3P + O_3P + M',
    fr={'Park', A=2.30e20, n=-1.00, T_a=128755.0, p_name='CO_X', p_mode='translation', s_p=0.5, q_name='CO_X', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=0.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=0.0,C_1D=0.0,C_1S=0.0,C_5So=0.0,C_plus=0.0,O_3P=0.0,O_1D=0.0,O_1S=0.0,O_plus=0.0,N_4So=0.0,N_2Do=0.0,N_2Po=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO_X', mode='vibration', model='TreanorMarrone', U=128755.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_a3 + M <=> C_3P + O_3P + M',
    fr={'Park', A=3.40e20, n=-1.00, T_a=58706.0, p_name='CO_a3', p_mode='translation', s_p=0.5, q_name='CO_a3', q_mode='vibration'},
    efficiencies={CO2=0.0,CO_X=0.0,CO_a3=0.0,CO_A=0.0,CO_plus=0.0,O2=0.0,N2_X=0.0,N2_A=0.0,N2_B=0.0,N2_C=0.0,NO=0.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=0.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=0.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO_a3', mode='vibration', model='TreanorMarrone', U=58706.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_a3 + M <=> C_3P + O_3P + M',
    fr={'Park', A=2.30e20, n=-1.00, T_a=58706.0, p_name='CO_a3', p_mode='translation', s_p=0.5, q_name='CO_a3', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=0.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=0.0,C_1D=0.0,C_1S=0.0,C_5So=0.0,C_plus=0.0,O_3P=0.0,O_1D=0.0,O_1S=0.0,O_plus=0.0,N_4So=0.0,N_2Do=0.0,N_2Po=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO_a3', mode='vibration', model='TreanorMarrone', U=58706.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_A + M <=> C_3P + O_3P + M',
    fr={'Park', A=3.40e20, n=-1.00, T_a=35126.0, p_name='CO_A', p_mode='translation', s_p=0.5, q_name='CO_A', q_mode='vibration'},
    efficiencies={CO2=0.0,CO_X=0.0,CO_a3=0.0,CO_A=0.0,CO_plus=0.0,O2=0.0,N2_X=0.0,N2_A=0.0,N2_B=0.0,N2_C=0.0,NO=0.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=0.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=0.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO_A', mode='vibration', model='TreanorMarrone', U=35126.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_A + M <=> C_3P + O_3P + M',
    fr={'Park', A=2.30e20, n=-1.00, T_a=35126.0, p_name='CO_A', p_mode='translation', s_p=0.5, q_name='CO_A', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=0.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=0.0,C_1D=0.0,C_1S=0.0,C_5So=0.0,C_plus=0.0,O_3P=0.0,O_1D=0.0,O_1S=0.0,O_plus=0.0,N_4So=0.0,N_2Do=0.0,N_2Po=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='CO_A', mode='vibration', model='TreanorMarrone', U=35126.0/3.0 } },
    ec={model='from thermo',iT=0}
}

-- 3.3 N2

reaction{
    'N2_X + M <=> N_4So + N_4So + M',
    fr={'Park', A=3.00e22, n=-1.60, T_a=113288.0, p_name='N2_X', p_mode='translation', s_p=0.5, q_name='N2_X', q_mode='vibration'},
    efficiencies={CO2=0.0,CO_X=0.0,CO_a3=0.0,CO_A=0.0,CO_plus=0.0,O2=0.0,N2_X=0.0,N2_A=0.0,N2_B=0.0,N2_C=0.0,NO=0.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=0.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=0.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='N2_X', mode='vibration', model='TreanorMarrone', U=113288.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_X + M <=> N_4So + N_4So + M',
    fr={'Park', A=7.01e21, n=-1.60, T_a=113288.0, p_name='N2_X', p_mode='translation', s_p=0.5, q_name='N2_X', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=0.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=0.0,C_1D=0.0,C_1S=0.0,C_5So=0.0,C_plus=0.0,O_3P=0.0,O_1D=0.0,O_1S=0.0,O_plus=0.0,N_4So=0.0,N_2Do=0.0,N_2Po=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='N2_X', mode='vibration', model='TreanorMarrone', U=113288.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_A + M <=> N_4So + N_4So + M',
    fr={'Park', A=3.00e22, n=-1.60, T_a=41057.0, p_name='N2_A', p_mode='translation', s_p=0.5, q_name='N2_A', q_mode='vibration'},
    efficiencies={CO2=0.0,CO_X=0.0,CO_a3=0.0,CO_A=0.0,CO_plus=0.0,O2=0.0,N2_X=0.0,N2_A=0.0,N2_B=0.0,N2_C=0.0,NO=0.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=0.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=0.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='N2_A', mode='vibration', model='TreanorMarrone', U=41057.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_A + M <=> N_4So + N_4So + M',
    fr={'Park', A=7.01e21, n=-1.60, T_a=41057.0, p_name='N2_A', p_mode='translation', s_p=0.5, q_name='N2_A', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=0.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=0.0,C_1D=0.0,C_1S=0.0,C_5So=0.0,C_plus=0.0,O_3P=0.0,O_1D=0.0,O_1S=0.0,O_plus=0.0,N_4So=0.0,N_2Do=0.0,N_2Po=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='N2_A', mode='vibration', model='TreanorMarrone', U=41057.0/3.0 } },
    ec={model='from thermo',iT=0}
}

-- note: seems strange that activation energy is higher for N2_B, need to check this
-- actually the dissociation potention is lower for the N2_B state

reaction{
    'N2_B + M <=> N_4So + N_4So + M',
    fr={'Park', A=3.00e22, n=-1.60, T_a=55175.0, p_name='N2_B', p_mode='translation', s_p=0.5, q_name='N2_B', q_mode='vibration'},
    efficiencies={CO2=0.0,CO_X=0.0,CO_a3=0.0,CO_A=0.0,CO_plus=0.0,O2=0.0,N2_X=0.0,N2_A=0.0,N2_B=0.0,N2_C=0.0,NO=0.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=0.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=0.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='N2_B', mode='vibration', model='TreanorMarrone', U=55175.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_B + M <=> N_4So + N_4So + M',
    fr={'Park', A=7.01e21, n=-1.60, T_a=55175.0, p_name='N2_B', p_mode='translation', s_p=0.5, q_name='N2_B', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=0.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=0.0,C_1D=0.0,C_1S=0.0,C_5So=0.0,C_plus=0.0,O_3P=0.0,O_1D=0.0,O_1S=0.0,O_plus=0.0,N_4So=0.0,N_2Do=0.0,N_2Po=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='N2_B', mode='vibration', model='TreanorMarrone', U=55175.0/3.0 } },
    ec={model='from thermo',iT=0}
}

-- 3.4 O2

reaction{
    'O2 + M <=> O_3P + O_3P + M',
    fr={'Park', A=1.0e22,  n=-1.5, T_a=59750.0, p_name='O2', p_mode='translation', s_p=0.5, q_name='O2', q_mode='vibration'},
    efficiencies={CO2=0.0,CO_X=0.0,CO_a3=0.0,CO_A=0.0,CO_plus=0.0,O2=0.0,N2_X=0.0,N2_A=0.0,N2_B=0.0,N2_C=0.0,NO=0.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=0.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=0.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='O2', mode='vibration', model='TreanorMarrone', U=59750.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'O2 + M <=> O_3P + O_3P + M',
    fr={'Park', A=2.0e21,  n=-1.5, T_a=59750.0, p_name='O2', p_mode='translation', s_p=0.5, q_name='O2', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=0.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=1.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=0.0,C_1D=0.0,C_1S=0.0,C_5So=0.0,C_plus=0.0,O_3P=0.0,O_1D=0.0,O_1S=0.0,O_plus=0.0,N_4So=0.0,N_2Do=0.0,N_2Po=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='O2', mode='vibration', model='TreanorMarrone', U=59750.0/3.0 } },
    ec={model='from thermo',iT=0}
}

-- not sure if we should set the efficiencies of N2_B and N2_C to zero in the above reaction due to these below

reaction{
    'O2 + N2_B <=> O_3P + O_3P + N2_X',
    fr={'Park', A=1.81e14,  n=0.0, T_a=0.0, p_name='O2', p_mode='translation', s_p=0.5, q_name='O2', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

reaction{
    'O2 + N2_C <=> O_3P + O_1S + N2_X',
    fr={'Park', A=1.81e14,  n=0.0, T_a=0.0, p_name='O2', p_mode='translation', s_p=0.5, q_name='O2', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

-- 3.5 NO

reaction{
    'NO + M <=> N_4So + O_3P + M',
    fr={'Park', A=9.63e14,  n=0.0, T_a=75210.0, p_name='NO', p_mode='translation', s_p=0.5, q_name='NO', q_mode='vibration'},
    efficiencies={CO2=1.0,CO_X=1.0,CO_a3=1.0,CO_A=1.0,CO_plus=0.0,O2=0.0,N2_X=0.0,N2_A=0.0,N2_B=0.0,N2_C=0.0,NO=1.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=1.0,C_1D=1.0,C_1S=1.0,C_5So=1.0,C_plus=0.0,O_3P=1.0,O_1D=1.0,O_1S=1.0,O_plus=0.0,N_4So=1.0,N_2Do=1.0,N_2Po=1.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='NO', mode='vibration', model='TreanorMarrone', U=75210.0/3.0 } },
    ec={model='from thermo',iT=0}
}

reaction{
    'NO + M <=> N_4So + O_3P + M',
    fr={'Park', A=1.45e15, n=0.0, T_a=75210.0, p_name='NO', p_mode='vibration', s_p=0.5, q_name='NO', q_mode='translation'},
    efficiencies={CO2=0.0,CO_X=0.0,CO_a3=0.0,CO_A=0.0,CO_plus=0.0,O2=1.0,N2_X=1.0,N2_A=1.0,N2_B=1.0,N2_C=1.0,NO=0.0,CN_X=0.0,CN_A=0.0,CN_B=0.0,C2_X=0.0,C2_d3=0.0,
                  C_3P=0.0,C_1D=0.0,C_1S=0.0,C_5So=0.0,C_plus=0.0,O_3P=0.0,O_1D=0.0,O_1S=0.0,O_plus=0.0,N_4So=0.0,N_2Do=0.0,N_2Po=0.0,N_plus=0.0,e_minus=0.0},
    chemistry_energy_coupling={ {species='NO', mode='vibration', model='TreanorMarrone', U=75210.0/3.0 } },
    ec={model='from thermo',iT=0}
}

-- 4. Electron-impact dissociation

-- 4.1 CN

reaction{
    'CN_X + e- <=> C_3P + N_4So + e-',
    fr={'Park', A=3.05e13, n=0.432, T_a=88966.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='CN_X', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_A + e- <=> C_3P + N_4So + e-',
    fr={'Park', A=1.04e14, n=0.464, T_a=75564.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='CN_A', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CN_B + e- <=> C_3P + N_4So + e-',
    fr={'Park', A=2.02e13, n=0.549, T_a=51576.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='CN_B', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

-- 4.2 CO

reaction{
    'CO_X + e- <=> C_3P + O_3P + e-',
    fr={'Park', A=8.65e12, n=0.367, T_a=129271.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='CO_X', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_a3 + e- <=> C_3P + O_3P + e-',
    fr={'Park', A=3.40e13, n=0.519, T_a=58742.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='CO_a3', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

reaction{
    'CO_A + e- <=> C_3P + O_3P + e-',
    fr={'Park', A=1.77e13, n=0.644, T_a=34864.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='CO_A', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

-- 4.3 N2

-- the A parameter for this reaction seems very low, but note that n is very large

reaction{
    'N2_X + e- <=> N_4So + N_4So + e-',
    fr={'Park', A=2.47e-9, n=6.16, T_a=113263.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='N2_X', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_A + e- <=> N_4So + N_4So + e-',
    fr={'Park', A=3.98e4, n=2.98, T_a=41670.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='N2_A', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_B + e- <=> N_4So + N_4So + e-',
    fr={'Park', A=2.71e0, n=3.73, T_a=55587.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='N2_B', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

reaction{
    'N2_C + e- <=> N_4So + N_4So + e-',
    fr={'Park', A=3.09e3, n=3.27, T_a=12893.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='N2_C', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

-- 4.4 O2

reaction{
    'O2 + e- <=> O_3P + O_3P + e-',
    fr={'Park', A=3.47e2,  n=3.52, T_a=59370.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='O2', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

-- 4.5 NO

reaction{
    'NO + e- <=> N_4So + O_3P + e-',
    fr={'Park', A=1.05e-2,  n=4.52, T_a=75390.0, p_name='e_minus', p_mode='translation', s_p=0.5, q_name='NO', q_mode='vibration'},
    ec={model='from thermo',iT=0}
}

-- 5. Neutral exchange

reaction{
   'N2_X + O_3P <=> NO + N_4So',
   fr={'Arrhenius', A=5.69e12,  n=0.42, T_a=42938.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O2 + N_4So <=> NO + O_3P',
   fr={'Arrhenius', A=2.49e12,  n=1.179, T_a=4005.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO_X + C_3P <=> C2_X + O_3P',
   fr={'Arrhenius', A=2.0e17,  n=-1.0, T_a=58000.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO_X + O_3P <=> O2 + C_3P',
   fr={'Arrhenius', A=3.9e13,  n=-0.18, T_a=69200.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO_X + N_4So <=> CN_X + O_3P',
   fr={'Arrhenius', A=1.0e14,  n=0.0, T_a=38600.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_X + C_3P <=> CN_X + N_4So',
   fr={'Arrhenius', A=5.24e13,  n=0.0, T_a=22600.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CN_X + O_3P <=> NO + C_3P',
   fr={'Arrhenius', A=1.6e13,  n=0.1, T_a=14600.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CN_X + C_3P <=> C2_X + N_4So',
   fr={'Arrhenius', A=5.0e13,  n=0.0, T_a=13000.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO_X + CO_X <=> C_3P + CO2',
   fr={'Arrhenius', A=2.3e9,  n=0.5, T_a=65710.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO2 + O_3P <=> O2 + CO_X',
   fr={'Arrhenius', A=2.1e13,  n=0.0, T_a=27800.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'C2_X + N2_X <=> CN_X + CN_X',
   fr={'Arrhenius', A=1.5e13,  n=0.0, T_a=21000.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO_X + N_4So <=> NO + C_3P',
   fr={'Arrhenius', A=2.9e11,  n=0.5, T_a=53630.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'NO + CO_X <=> CO2 + N_4So',
   fr={'Arrhenius', A=4.6e8,  n=0.5, T_a=12070.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'NO + NO <=> N2_X + O2',
   fr={'Arrhenius', A=3.07e11,  n=0.0, T_a=33660.0},
   ec={model='from thermo',iT=0}
}

-- 6. Neutral exchange involving excited states

reaction{
   'CO_a3 + O_3P <=> O2 + C_3P',
   fr={'Arrhenius', A=4.0e13,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_A + CO_X <=> N2_X + CO_a3',
   fr={'Arrhenius', A=1.2e13,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_A + CO_X <=> N2_X + CO_A',
   fr={'Arrhenius', A=1.1e10,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_A + O_3P <=> NO + N_2Do',
   fr={'Arrhenius', A=4.21e12,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N_2Do + O2 <=> NO + O_3P',
   fr={'Arrhenius', A=3.13e12,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N_2Do + NO <=> N2_X + O_3P',
   fr={'Arrhenius', A=1.08e14,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N_2Po + O2 <=> NO + O_3P',
   fr={'Arrhenius', A=1.08e14,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N_2Po + NO <=> N2_X + O_3P',
   fr={'Arrhenius', A=1.81e13,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O_1D + NO <=> O2 + N_4So',
   fr={'Arrhenius', A=1.02e14,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

-- 7. Quenching

reaction{
   'N2_A + O_3P <=> N2_X + O_1S',
   fr={'Arrhenius', A=1.26e13,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_A + N_4So <=> N2_X + N_4So',
   fr={'Arrhenius', A=1.20e12,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_A + N_4So <=> N2_X + N_2Po',
   fr={'Arrhenius', A=1.08e15,  n=-0.67, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_A + N2_X <=> N2_X + N2_X',
   fr={'Arrhenius', A=1.81e8,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_A + N2_A <=> N2_X + N2_B',
   fr={'Arrhenius', A=1.81e14,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_A + N2_A <=> N2_X + N2_C',
   fr={'Arrhenius', A=9.03e13,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_B + N2_X <=> N2_X + N2_X',
   fr={'Arrhenius', A=1.20e12,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N2_B + NO <=> N2_A + NO',
   fr={'Arrhenius', A=1.44e14,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N_2Do + O_3P <=> N_4So + O_1D',
   fr={'Arrhenius', A=2.41e11,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N_2Do + N2_X <=> N_4So + N2_X',
   fr={'Arrhenius', A=3.61e9,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N_2Po + N_4So <=> N_2Do + N_4So',
   fr={'Arrhenius', A=1.08e12,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N_2Po + N2_X <=> N_4So + N2_X',
   fr={'Arrhenius', A=1.20e6,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O_1D + O_3P <=> O_3P + O_3P',
   fr={'Arrhenius', A=4.82e12,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O_1D + O2 <=> O_3P + O2',
   fr={'Arrhenius', A=3.85e12,  n=0.0, T_a=67.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O_1D + N2_X <=> O_3P + N2_X',
   fr={'Arrhenius', A=1.38e12,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O_1S + O_3P <=> O_1D + O_1D',
   fr={'Arrhenius', A=3.01e13,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O_1S + N_4So <=> O_3P + N_4So',
   fr={'Arrhenius', A=6.02e11,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O_1S + O2 <=> O_3P + O2',
   fr={'Arrhenius', A=7.83e11,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O_1S + N2_X <=> O_3P + N2_X',
   fr={'Arrhenius', A=6.02e6,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O_1S + NO <=> O_3P + NO',
   fr={'Arrhenius', A=1.75e14,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O_1S + NO <=> O_1D + NO',
   fr={'Arrhenius', A=1.75e14,  n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

-- 8. Associative ionisation

reaction{
   'CO+ + e- <=> C_3P + O_3P',
   fr={'Park', A=1.86e18, n=-0.48, T_a=0.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   ec={model='from thermo',iT=0}
}

reaction{
   'N_4So + O_3P <=> NO+ + e-',
   fr={'Arrhenius', A=8.8e8, n=1.0, T_a=31900.0},
   ec={model='from thermo',iT=0}
   -- New reaction: from Park (1994) scheme
}

reaction{
   'O_3P + O_3P <=> O2+ + e-',
   fr={'Arrhenius', A=7.1e2, n=2.7, T_a=80600.0},
   ec={model='from thermo',iT=0}
   -- New reaction: from Park (1994) scheme
}


-- 9. Electron impact ionisation 

reaction{
   'C_3P + e- <=> C+ + e- + e-',
   fr={'Park', A=5.85e27, n=-2.074, T_a=127510.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=127510.0} },
   ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

reaction{
   'N_4So + e- <=> N+ + e- + e-',
   fr={'Park', A=1.93e31, n=-2.856, T_a=168970.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=168970.0} },
   ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

-- NOTE: the A parameter in the following reaction has been corrected from 8.25e37 to 8.25e27, based on comparison with the park rate

reaction{
   'O_3P + e- <=> O+ + e- + e-',
   fr={'Park', A=8.25e27, n=-2.237, T_a=157840.0, p_name='e_minus', p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   chemistry_energy_coupling={ {species='e_minus', mode='translation', model='electron impact ionization', T_I=157840.0} },
   ec={model='from thermo',iT=-1,species='e_minus',mode='translation'}
}

-- 10. Charge exchange reactions

reaction{
   'CO_X + C+ <=> CO+ + C_3P',
   fr={'Arrhenius', A=1.0e13, n=0.0, T_a=31400.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'CO+ + O_3P <=> CO_X + O+',
   fr={'Arrhenius', A=8.43e13, n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N+ + O_3P <=> N_4So + O+',
   fr={'Arrhenius', A=6.02e11, n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'O+ + N_2Do <=> N+ + O_3P',
   fr={'Arrhenius', A=7.83e13, n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'NO+ + C_3P <=> NO + C+',
   fr={'Arrhenius', A=1.0e13,  n=0.0, T_a=23200.0},
   ec={model='from thermo',iT=0}
   -- New reaction: from Park (1994) scheme
}

reaction{
   'O2+ + O_3P <=> O+ + O2',
   fr={'Arrhenius', A=4.0e12,  n=-0.09, T_a=18000.0},
   ec={model='from thermo',iT=0}
   -- New reaction: from Park (1994) scheme
}

reaction{
   'NO+ + N_4So <=> O+ + N2_X',
   fr={'Arrhenius', A=3.4e13,  n=-1.08, T_a=12800.0},
   ec={model='from thermo',iT=0}
   -- New reaction: from Park (1994) scheme
}

reaction{
   'NO+ + O_3P <=> O2+ + N_4So',
   fr={'Arrhenius', A=7.2e12,  n=0.29, T_a=48600.0},
   ec={model='from thermo',iT=0}
   -- New reaction: from Park (1994) scheme
}

reaction{
   'O2 + C+ <=> O2+ + C_3P',
   fr={'Arrhenius', A=1.0e13,  n=0.0, T_a=9400.0},
   ec={model='from thermo',iT=0}
   -- New reaction: from Park (1994) scheme
}

-- 11. Ion exchange reactions

reaction{
   'N+ + O2 <=> O+ + NO',
   fr={'Arrhenius', A=1.69e13, n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

reaction{
   'N+ + NO <=> N2_X + O+',
   fr={'Arrhenius', A=6.02e11, n=0.0, T_a=0.0},
   ec={model='from thermo',iT=0}
}

-- 12. Radiative transitions

-- 12.1 C2

reaction{
   'C2_d3 <=> C2_X',
   fr={'Arrhenius', A=9.3e6, n=0.0, T_a=0.0},
   br={'Arrhenius', A=1.0e-100, n=0.0, T_a=0.0},
   -- Zalogin
}


-- 12.1 CN

reaction{
   'CN_B <=> CN_A',
   fr={'Arrhenius', A=9.3e5, n=0.0, T_a=0.0},
   br={'Arrhenius', A=1.0e-100, n=0.0, T_a=0.0},
   -- From my thesis
}

reaction{
   'CN_A <=> CN_X',
   fr={'Arrhenius', A=1.7e6, n=0.0, T_a=0.0},
   br={'Arrhenius', A=1.0e-100, n=0.0, T_a=0.0},
   -- From my thesis
}

reaction{
   'CN_B <=> CN_X',
   fr={'Arrhenius', A=2.0e8, n=0.0, T_a=0.0},
   br={'Arrhenius', A=1.0e-100, n=0.0, T_a=0.0},
   -- From my thesis
}

-- 12.2 CO

-- this transition is optically thick
-- reaction{
--    'CO_A <=> CO_X',
--    fr={'Arrhenius', A=1.4e9, n=0.0, T_a=0.0},
--    br={'Arrhenius', A=1.0e-100, n=0.0, T_a=0.0},
   -- From my thesis
-- }

-- 12.3 N2

reaction{
   'N2_B <=> N2_A',
   fr={'Arrhenius', A=1.4e5, n=0.0, T_a=0.0},
   br={'Arrhenius', A=1.0e-100, n=0.0, T_a=0.0},
   -- Chernyi via CJ
}

reaction{
   'N2_C <=> N2_B',
   fr={'Arrhenius', A=2.6e7, n=0.0, T_a=0.0},
   br={'Arrhenius', A=1.0e-100, n=0.0, T_a=0.0},
   -- Pancheshnyi via CJ
}
