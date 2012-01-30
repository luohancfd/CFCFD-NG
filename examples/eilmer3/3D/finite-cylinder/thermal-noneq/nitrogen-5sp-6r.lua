-- Author: Daniel F. Potter
-- Date: 29-Nov-2011
-- Place: DLR Göttingen, Germany
--
-- Two-temperature version of the ionising nitrogen
-- reaction scheme.
-- 
-- Updated from the work by RJG and DFP as found in
-- lib/gas_models2/input_files/nitrogen/nitrogen-5sp-6r.py
--
-- Note: Based on Dan's comments, I've only included
-- the Goekcen rates at present.
--
-- Reference:
-- Goekcen (2004)
-- N2-CH4-Ar Chemical Kinetic Model for Simulations of 
-- Atmospheric Entry to Titan
-- AIAA Paper 2004-2469
--

scheme_t = {
    update = "chemical kinetic ODE MC",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

reaction{
   'N2 + N2 <=> N + N + N2',
   fr={'Park', A=7.0e21, n=-1.6, T_a=113200.0, p_name='N2', p_mode='vibration', 
       s_p=0.3, q_name='N2', q_mode='translation'},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + N <=> N + N + N',
   fr={'Park', A=3.0e22, n=-1.6, T_a=113200.0, p_name='N2', p_mode='vibration', 
       s_p=0.3, q_name='N', q_mode='translation'},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + e- <=> N + N + e-',
   fr={'Park', A=3.0e24, n=-1.6, T_a=113200.0, p_name='N2', p_mode='vibration', 
       s_p=0.3, q_name='e_minus', q_mode='translation'},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N + N <=> N2+ + e-',
   fr={'Arrhenius', A=4.40e7, n=1.5, T_a=67500.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + N+ <=> N2+ + N',
   fr={'Arrhenius', A=1.0e12, n=0.5, T_a=12200.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N + e- <=> N+ + e- + e-',
   fr={'Park', A=2.50e34, n=-3.82, T_a=168600.0, p_name='e_minus', 
       p_mode='translation', s_p=1.0, q_name='NA', q_mode='NA'},
   ec={model='from CEA curves', iT=-1, species='e_minus', mode='translation'}
}

