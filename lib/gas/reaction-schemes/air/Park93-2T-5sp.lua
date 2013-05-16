-- Park93-2T-5sp.lua
--
-- The single temperature, neutral reactions from:
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
--   16-May-2012: - Corrected the third-body efficiencies in NO dissociation reaction
--
-- ODE SOLVER
-- Chemical kinetic ODE MC

scheme_t = {
    update = "chemical kinetic ODE MC",
    temperature_limits = {
        lower = 290.0,
        upper = 100000.0
    },
    error_tolerance = 0.000000001
}

-- Dissociation reactions

reaction{
   'N2 + M <=> N + N + M',
   fr={'Park', A=7.0e21, n=-1.60, T_a=113200.0, p_name='N2', p_mode='vibration', s_p=0.5, q_name='N2', q_mode='translation'},
   efficiencies={N2=1.0,O2=1.0,NO=1.0,N=4.2857,O=4.2857},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O2 + M <=> O + O + M',
   fr={'Park', A=2.0e21, n=-1.5, T_a=59500.0, p_name='O2', p_mode='vibration', s_p=0.5, q_name='O2', q_mode='translation'},
   efficiencies={N2=1.0,O2=1.0,NO=1.0,N=5.0,O=5.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'NO + M <=> N + O + M',
   fr={'Park', A=5.0e15, n=0.0, T_a=75500.0, p_name='NO', p_mode='vibration', s_p=0.5, q_name='NO', q_mode='translation'},
   efficiencies={N2=1.0,O2=1.0,NO=1.0,N=22.0,O=22.0},
   ec={model='from CEA curves',iT=0}
}

-- Neutral exchange reactions

reaction{
   'NO + O <=> O2 + N',
   fr={'Arrhenius', A=8.4e12, n=0.0, T_a=19450.0},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'N2 + O <=> NO + N',
   fr={'Arrhenius', A=6.4e17, n=-1.0, T_a=38400.0},
   ec={model='from CEA curves',iT=0}
}

