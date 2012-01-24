-- Jachimowski_88.py
--
-- This file provides a chemical kinetic description
-- of hydrogen combustion in air.
--
-- The numbering of reactions in this file corresponds
-- to Table 1 in Jachimowski (1998)
--
-- Reference:
--  Jachimowksi C.J. (1998)
--  An Analytical Study of the Hydrogen-Air Reaction
--  Mechanism wih Application to Scramjet Combustion
--
--  Modified model...
--  Oldenborg, R., Chinitz, W., Friedman, M., Jaffe, R.,
--  Jachimowski, C., Rabinowitz, M. and Schott, G. (1990)
--  Hypersonic Combustion Kinetics: Status report of the
--  rate constant committee, NASP High Speed Propulsion
--  Technology Team.
--  NASP Rate Constant Committee, NASP TM-1107
--
-- This file prepared by..
-- Rowan J Gollan
-- 15-May-2006
--
-- Adapted from python file Jachimowski_88.py
-- Fabian Zander
-- 19-Oct-2010
--
-- Species included: H, H2, O, O2, N, N2, OH, NO, H2O, HO2, NO2, H2O2, HNO

MODIFIED = false
O2_only = false

if MODIFIED then

	reaction{
	   'H2 + O2 <=> H + HO2',
	   fr={'Arrhenius', A=1.00e14, n=0.0, T_a=56000.0/1.987},
	   label='r1'
	}
else
	reaction{
	   'H2 + O2 <=> OH + OH',
	   fr={'Arrhenius', A=1.70e13, n=0.0, T_a=48000.0/1.987},
	   label='r1'
	}
end

reaction{
   'H + O2 <=> OH + O',
   fr={'Arrhenius', A=2.6e14, n=0.0, T_a=16800.0/1.987},
   label='r2'
}

reaction{
   'O + H2 <=> OH + H',
   fr={'Arrhenius', A=1.80e10, n=1.0, T_a=8900.0/1.987},
   label='r3'
}

reaction{
   'OH + H2 <=> H2O + H',
   fr={'Arrhenius', A=2.20e13, n=0.0, T_a=5150.0/1.987},
   label='r4'
}

reaction{
   'OH + OH <=> H2O + O',
   fr={'Arrhenius', A=6.30e12, n=0.0, T_a=1090.0/1.987},
   label='r5'
}

reaction{
   'H + OH + M <=> H2O + M',
   fr={'Arrhenius', A=2.20e22, n=-2.0, T_a=0.0/1.987},
   label='r6',
   efficiencies={['H2O']=6.0}
}

reaction{
   'H + H + M <=> H2 + M',
   fr={'Arrhenius', A=6.40e17, n=-1.0, T_a=0.0/1.987},
   label='r7',
   efficiencies={['H2']=2.0,['H2O']=6.0}
}

reaction{
   'H + O + M <=> OH + M',
   fr={'Arrhenius', A=6.00e16, n=-0.60, T_a=0.0/1.987},
   label='r8',
   efficiencies={['H2O']=5.0}
}

reaction{
   'H + O2 + M <=> HO2 + M',
   fr={'Arrhenius', A=2.10e15, n=0.00, T_a=-1000.0/1.987},
   label='r9',
   efficiencies={['H2']=2.0, ['H2O']=16.0}
}

reaction{
   'HO2 + H <=> H2 + O2',
   fr={'Arrhenius', A=1.30e13, n=0.00, T_a=0.0/1.987},
   label='r10',
}

reaction{
   'HO2 + H <=> OH + OH',
   fr={'Arrhenius', A=1.40e14, n=0.00, T_a=1080.0/1.987},
   label='r11',
}

reaction{
   'HO2 + H <=> H2O + O',
   fr={'Arrhenius', A=1.0e13, n=0.00, T_a=1080.0/1.987},
   label='r12',
}

reaction{
   'HO2 + O <=> O2 + OH',
   fr={'Arrhenius', A=1.50e13, n=0.00, T_a=950.0/1.987},
   label='r13',
}

reaction{
   'HO2 + OH <=> H2O + O2',
   fr={'Arrhenius', A=8.00e12, n=0.00, T_a=0.0/1.987},
   label='r14',
}

reaction{
   'HO2 + HO2 <=> H2O2 + O2',
   fr={'Arrhenius', A=2.00e12, n=0.00, T_a=0.0/1.987},
   label='r15',
}

reaction{
   'H + H2O2 <=> H2 + HO2',
   fr={'Arrhenius', A=1.40e12, n=0.00, T_a=3600.0/1.987},
   label='r16',
}

reaction{
   'O + H2O2 <=> OH + HO2',
   fr={'Arrhenius', A=1.40e13, n=0.00, T_a=6400.0/1.987},
   label='r17',
}

reaction{
   'OH + H2O2 <=> H2O + HO2',
   fr={'Arrhenius', A=6.10e12, n=0.00, T_a=1430.0/1.987},
   label='r18',
}

reaction{
   'M + H2O2 <=> OH + OH + M',
   fr={'Arrhenius', A=1.20e17, n=0.00, T_a=45500.0/1.987},
   efficiencies={['H2O']=15.0},
   label='r19',
}

reaction{
   'O + O + M <=> O2 + M',
   fr={'Arrhenius', A=6.00e17, n=0.00, T_a=-1800.0/1.987},
   label='r20',
}


reaction{
	'N + N + M <=> N2 + M',
	fr={'Arrhenius', A=2.80e17, n=-0.75, T_a=0.0/1.987},
	label='r21',
}

reaction{
	'N + O2 <=> NO + O',
	fr={'Arrhenius', A=6.40e9, n=1.00, T_a=6300.0/1.987},
	label='r22',
}

reaction{
	'N + NO <=> N2 + O',
	fr={'Arrhenius', A=1.60e13, n=0.00, T_a=0.0/1.987},
	label='r23',
}

reaction{
	'N + OH <=> NO + H',
	fr={'Arrhenius', A=6.30e11, n=0.50, T_a=0.0/1.987},
	label='r24',
}

reaction{
	'H + NO + M <=> HNO + M',
	fr={'Arrhenius', A=5.40e15, n=0.00, T_a=-600.0/1.987},
	label='r25',
}

reaction{
	'H + HNO <=> NO + H2',
	fr={'Arrhenius', A=4.80e12, n=0.00, T_a=0.0/1.987},
	label='r26',
}

reaction{
	'O + HNO <=> NO + OH',
	fr={'Arrhenius', A=5.00e11, n=0.50, T_a=0.0/1.987},
	label='r27',
}

reaction{
	'OH + HNO <=> NO + H2O',
	fr={'Arrhenius', A=3.60e13, n=0.00, T_a=0.0/1.987},
	label='r28',
}

reaction{
	'HO2 + HNO <=> NO + H2O2',
	fr={'Arrhenius', A=2.00e12, n=0.00, T_a=0.0/1.987},
	label='r29',
}

reaction{
	'HO2 + NO <=> NO2 + OH',
	fr={'Arrhenius', A=3.40e12, n=0.00, T_a=-260.0/1.987},
	label='r30',
}

reaction{
	'H + NO2 <=> NO + OH',
	fr={'Arrhenius', A=3.50e14, n=0.00, T_a=1500.0/1.987},
	label='r31',
}

reaction{
	'O + NO2 <=> NO + O2',
	fr={'Arrhenius', A=1.00e13, n=0.00, T_a=600.0/1.987},
	label='r32',
}

reaction{
	'M + NO2 <=> NO + O + M',
	fr={'Arrhenius', A=1.16e16, n=0.00, T_a=66000.0/1.987},
	label='r33',
}

if O2_only then
	select_reactions_by_label({'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18', 'r19', 'r20'})
end
