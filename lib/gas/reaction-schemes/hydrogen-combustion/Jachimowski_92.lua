-- Jachimowski_92.lua
--
-- This file provides a chemical kinetic description
-- of hydrogen combustion in air.
--
-- The numbering of reactions in this file corresponds
-- to Table II in Jachimowski (1992)
--
-- Reference:
--  Jachimowksi C.J. (1992)
--  An Analysis of Combustion Studies in Shock Expansion
--  Tunnels and Reflected Shock Tunnels
--  NASA TP-3224
--
-- This file prepared by
-- Fabian Zander
-- 29-Oct-2010
--
-- Species included: H, H2, O, O2, N, N2, OH, NO, H2O, HO2, NO2, H2O2, HNO

reaction{
	'H2 + O2 <=> HO2 + H',
	fr={'Arrhenius', A=7.00e13, n=0.0, T_a=56800.0/1.987},
	label='r1'
}

reaction{
	'H + O2 <=> OH + O',
	fr={'Arrhenius', A=2.2e14, n=0.0, T_a=16800.0/1.987},
	label='r2'
}

reaction{
	'O + H2 <=> OH + H',
	fr={'Arrhenius', A=5.06e4, n=2.67, T_a=6290.0/1.987},
	label='r3'
}

reaction{
	'OH + H2 <=> H2O + H',
	fr={'Arrhenius', A=2.16e8, n=1.51, T_a=3430.0/1.987},
	label='r4'
}

reaction{
	'OH + OH <=> H2O + O',
	fr={'Arrhenius', A=1.50e9, n=1.14, T_a=0.0/1.987},
	label='r5'
}

reaction{
	'H + OH + M <=> H2O + M',
	fr={'Arrhenius', A=8.62e21, n=-2.0, T_a=0.0/1.987},
	label='r6',
	efficiencies={['H2O']=16.0, ['H2']=2.5, ['HO2']=4.0}
}

reaction{
	'H + H + M <=> H2 + M',
	fr={'Arrhenius', A=7.30e17, n=-1.0, T_a=0.0/1.987},
	label='r7',
	efficiencies={['H2O']=16.0, ['H2']=2.5, ['HO2']=4.0}
}

reaction{
	'H + O + M <=> OH + M',
	fr={'Arrhenius', A=2.60e16, n=-0.60, T_a=0.0/1.987},
	label='r8',
	efficiencies={['H2O']=16.0, ['H2']=2.5, ['HO2']=4.0}
}

reaction{
	'O + O + M <=> O2 + M',
	fr={'Arrhenius', A=1.10e17, n=-1.00, T_a=0.0/1.987},
	label='r9',
	efficiencies={['H2O']=16.0, ['H2']=2.5, ['HO2']=4.0}
}

reaction{
	'H + O2 + M <=> HO2 + M',
	fr={'Arrhenius', A=2.30e18, n=-1.00, T_a=0.0/1.987},
	label='r10',
	efficiencies={['H2O']=16.0, ['H2']=2.5, ['HO2']=4.0}
}

reaction{
	'HO2 + H <=> OH + OH',
	fr={'Arrhenius', A=1.50e14, n=0.00, T_a=1000.0/1.987},
	label='r11',
}

reaction{
	'HO2 + O <=> O2 + OH',
	fr={'Arrhenius', A=2.00e13, n=0.00, T_a=0.0/1.987},
	label='r12',
}

reaction{
	'HO2 + OH <=> H2O + O2',
	fr={'Arrhenius', A=2.0e13, n=0.00, T_a=0.0/1.987},
	label='r13',
}

reaction{
	'HO2 + HO2 <=> H2O2 + O2',
	fr={'Arrhenius', A=2.00e12, n=0.00, T_a=0.0/1.987},
	label='r14',
}

reaction{
	'H + H2O2 <=> H2 + HO2',
	fr={'Arrhenius', A=1.70e12, n=0.00, T_a=3780.0/1.987},
	label='r15',
}

reaction{
	'H + H2O2 <=> OH + H2O',
	fr={'Arrhenius', A=1.00e13, n=0.00, T_a=3580.0/1.987},
	label='r16',
}

reaction{
	'O + H2O2 <=> OH + HO2',
	fr={'Arrhenius', A=2.80e13, n=0.00, T_a=6400.0/1.987},
	label='r17',
}

reaction{
	'OH + H2O2 <=> H2O + HO2',
	fr={'Arrhenius', A=7.00e12, n=0.00, T_a=1435.0/1.987},
	label='r18',
}

reaction{
	'OH + OH + M <=> H2O2 + M',
	fr={'Arrhenius', A=1.60e22, n=-2.00, T_a=0.0/1.987},
	label='r19',
	efficiencies={['H2O']=16.0, ['H2']=2.5, ['HO2']=4.0}
}

reaction{
	'N + N + M <=> N2 + M',
	fr={'Arrhenius', A=2.80e17, n=-0.8, T_a=0.0/1.987},
	label='r20',
	efficiencies={['H2O']=16.0, ['H2']=2.5, ['HO2']=4.0}
}

reaction{
	'N + O2 <=> NO + O',
	fr={'Arrhenius', A=6.40e9, n=1.00, T_a=6300.0/1.987},
	label='r21',
}

reaction{
	'N + NO <=> N2 + O',
	fr={'Arrhenius', A=1.60e13, n=0.00, T_a=0.0/1.987},
	label='r22',
}

reaction{
	'N + OH <=> NO + H',
	fr={'Arrhenius', A=6.30e11, n=0.50, T_a=0.0/1.987},
	label='r23',
}

reaction{
	'H + NO + M <=> HNO + M',
	fr={'Arrhenius', A=5.40e15, n=0.00, T_a=-600.0/1.987},
	label='r24',
	efficiencies={['H2O']=16.0, ['H2']=2.5, ['HO2']=4.0}
}

reaction{
	'H + HNO <=> NO + H2',
	fr={'Arrhenius', A=4.80e12, n=0.00, T_a=0.0/1.987},
	label='r25',
}

reaction{
	'O + HNO <=> NO + OH',
	fr={'Arrhenius', A=5.00e11, n=0.50, T_a=0.0/1.987},
	label='r26',
}

reaction{
	'OH + HNO <=> NO + H2O',
	fr={'Arrhenius', A=3.60e13, n=0.00, T_a=0.0/1.987},
	label='r27',
}

reaction{
	'HO2 + HNO <=> NO + H2O2',
	fr={'Arrhenius', A=2.00e12, n=0.00, T_a=0.0/1.987},
	label='r28',
}

reaction{
	'HO2 + NO <=> NO2 + OH',
	fr={'Arrhenius', A=3.40e12, n=0.00, T_a=-260.0/1.987},
	label='r29',
}

reaction{
	'HO2 + NO <=> HNO + O2',
	fr={'Arrhenius', A=2.00e11, n=0.00, T_a=1000.0/1.987},
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
	efficiencies={['H2O']=16.0, ['H2']=2.5, ['HO2']=4.0}
}
