-- Bittker_Scullin.py
--
-- This file provides a chemical kinetic description
-- of hydrogen combustion in air.
--
-- The numbering of reactions in this file corresponds
-- to the scheme on p.85 of Bittker and Scullin
--
-- Reference:
--  D. A. Bittker and V. J. Scullin
--  General Chemical Kinetics Computer Program for Static
--  and Flow Reactions, with Application to Combustion
--  and Shock Tube Kinetics
--  NASA TN D-6586 1972
--
-- This file prepared by..
-- Fabian Zander
-- 28th Oct 2010
--
-- a 9-species, 18-reaction scheme
-- Species included: H, H2, O, O2, N2, OH, H2O, HO2, H2O2

reaction{
   'H2 + OH <=> H2O + H',
   fr={'Arrhenius', A=2.10e13, n=0.0, T_a=5100.0/1.987},
   label='r1'
}

reaction{
   'H + O2 <=> OH + O',
   fr={'Arrhenius', A=1.25e14, n=0.0, T_a=16300.0/1.987},
   label='r2'
}

reaction{
   'O + H2 <=> OH + H',
   fr={'Arrhenius', A=2.95e13, n=1.0, T_a=9800.0/1.987},
   label='r3'
}

reaction{
   'H + O2 + M <=> HO2 + M',
   fr={'Arrhenius', A=1.59e15, n=0.00, T_a=-1000.0/1.987},
   label='r4',
   efficiencies={['H2']=5.0, ['H2O']=32.5, ['O2']=2.0, ['N2']=2.0}
}

reaction{
   'H + H + M <=> H2 + M',
   fr={'Arrhenius', A=1.0e18, n=-1.0, T_a=0.0/1.987},
   label='r5',
   efficiencies={['H2']=5.0,['H2O']=15.0, ['O2']=2.0, ['N2']=2.0}
}

reaction{
   'H2 + HO2 <=> H2O2 + H',
   fr={'Arrhenius', A=9.6e12, n=0.00, T_a=24000.0/1.987},
   label='r6',
}

reaction{
   'M + H2O2 <=> OH + OH + M',
   fr={'Arrhenius', A=1.17e17, n=0.00, T_a=45500.0/1.987},
   efficiencies={['H2O2']=6.6, ['H2']=2.3, ['H2O']=6.0, ['O2']=0.78},
   label='r7',
}

reaction{
   'HO2 + H <=> OH + OH',
   fr={'Arrhenius', A=7.0e13, n=0.00, T_a=0.0/1.987},
   label='r8',
}

reaction{
   'H + OH + M <=> H2O + M',
   fr={'Arrhenius', A=7.50e23, n=-2.6, T_a=0.0/1.987},
   label='r9',
   efficiencies={['H2']=4.0, ['H2O']=20.0, ['O2']=1.6, ['N2']=1.6}
}

reaction{
   'O + O + M <=> O2 + M',
   fr={'Arrhenius', A=1.38e18, n=-1.00, T_a=340.0/1.987},
   label='r10',
}

reaction{
   'O + H2O <=> OH + OH',
   fr={'Arrhenius', A=5.75e13, n=0.00, T_a=18000.0/1.987},
   label='r11',
}

reaction{
   'H2 + O2 <=> OH + OH',
   fr={'Arrhenius', A=1.00e13, n=0.0, T_a=43000.0/1.987},
   label='r12'
}

reaction{
   'HO2 + OH <=> H2O + O2',
   fr={'Arrhenius', A=6.30e12, n=0.00, T_a=0.0/1.987},
   label='r13',
}

reaction{
   'HO2 + O <=> O2 + OH',
   fr={'Arrhenius', A=6.00e12, n=0.00, T_a=0.0/1.987},
   label='r14',
}

reaction{
   'HO2 + HO2 <=> H2O2 + O2',
   fr={'Arrhenius', A=1.80e12, n=0.00, T_a=0.0/1.987},
   label='r15',
}

reaction{
   'OH + H2O2 <=> H2O + HO2',
   fr={'Arrhenius', A=1.00e13, n=0.00, T_a=1800.0/1.987},
   label='r16',
}

reaction{
   'O + H2O2 <=> OH + HO2',
   fr={'Arrhenius', A=8.00e13, n=0.00, T_a=1000.0/1.987},
   label='r17',
}

reaction{
   'H + H2O2 <=> H2O + OH',
   fr={'Arrhenius', A=3.18e14, n=0.00, T_a=9000.0/1.987},
   label='r18',
}
