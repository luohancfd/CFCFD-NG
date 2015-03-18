-- Author: Rowan J. Gollan
-- Date: 2015-03-18
--
-- Reference:
-- Burke, M.P., Chaos, M., Ju, Y., Dryer, F.L. and Klippenstein, S.J. (2012)
-- Comprehensive H2/O2 kinetic model for high-pressure combustion,
-- International Journal of Chemical Kinetics, 44:7, pp.444-474
--
-- Note:
-- I have added notes comparing the rate coefficients to those
-- listed in the Stanford mechanism. These mechansisms are
-- very similar in many ways. When I note something as "identical"
-- I mean to within significant figures (so you might see some
-- small differences in the numerical value).


S = 1.0/1.987

reaction{
   'H + O2 <=> O + OH',
   fr={'Arrhenius', A=1.04e14, n=0.0, T_a=1.531e4*S},
   label='r1',
   note='very similar to Stanford'
}

reaction{
   'O + H2 <=> H + OH',
   fr={'Arrhenius', A=3.82e12, n=0.0, T_a=7.948e3*S},
   label='r2a',
   note='identical to Stanford r9a'
}
reaction{
   'O + H2 <=> H + OH',
   fr={'Arrhenius', A=8.79e12, n=0.0, T_a=1.917e4*S},
   label='r2b',
   note='identical to Stanford r9b'
}

reaction{
   'H2 + OH <=> H2O + H',
   fr={'Arrhenius', A=2.16e8, n=1.51, T_a=3.430e3*S},
   label='r3',
   note='very similar to Stanford, r10'
}

reaction{
   'OH + OH <=> O + H2O',
   fr={'Arrhenius', A=3.34e4, n=2.42, T_a=-1.930e3*S},
   label='r4',
   note='similar to Stanford, r8'
}

reaction{
   'H2 + M <=> H + H + M',
   fr={'Arrhenius', A=4.58e19, n=-1.40, T_a=1.040e5*S},
   efficiencies={H2=2.5,H2O=12.0},
   label='r5',
   note='similar to Stanford (r18) when applied to O2, differs w.r.t H2'
}

reaction{
   'O + O + M <=> O2 + M',
   fr={'Arrhenius', A=6.16e15, n=-0.50, T_a=0.0},
   efficiencies={H2=2.5,H2O=12.0},
   label='r6',
   note='identical to Stanford, r19 (except that Stanford chem.inp give H2O2 efficiency'
}

reaction{
   'O + H + M <=> OH + M',
   fr={'Arrhenius', A=4.71e18, n=-1.0, T_a=0.0},
   efficiencies={H2=2.5,H2O=12.0},
   label='r7',
   note='identical to Stanford, r20 (except that Stanford chem.inp give H2O2 efficiency'
}

reaction{
   'H2O + M <=> H + OH + M',
   fr={'Arrhenius', A=6.06e27, n=-3.32, T_a=1.208e5*S},
   efficiencies={H2=3.0,H2O=0.0},
   label='r8a',
   note='identical to Stanford, r7a (except for O2 efficiency)'
}
reaction{
   'H2O + H2O <=> H + OH + H2O',
   fr={'Arrhenius', A=1.01e26, n=-2.44, T_a=1.202e5*S},
   label='r8b',
   note='identical to Stanford, r7b'
}

reaction{
   'H + O2 (+ M) <=> HO2 (+ M)',
   fr={'pressure dependent',
       k_inf={A=4.65e12, n=0.44, T_a=0.0},
       k_0={A=6.37e20, n=-1.72, T_a=5.250e2*S},
       Troe={F_cent=0.5}
   },
   efficiencies={H2=2.0,H2O=14.0},
   label='r9',
   note=[[The first option is selected for N2 as bath gas.
          This differs somewhat from Stanford r2]]
}

reaction{
   'HO2 + H <=> H2 + O2',
   fr={'Arrhenius', A=2.75e6, n=2.09, T_a=-1.451e3*S},
   label='r10',
   note='differs from Stanford, r13'
}

reaction{
   'HO2 + H <=> OH + OH',
   fr={'Arrhenius', A=7.08e13, n=0.0, T_a=2.95e2*S},
   label='r11',
   note='identical to Stanford, r11'
}

reaction{
   'HO2 + O <=> O2 + OH',
   fr={'Arrhenius', A=2.85e10, n=1.0, T_a=-7.239e2*S},
   label='r12',
   note='differs from Stanford, r14'
}

reaction{
   'HO2 + OH <=> H2O + O2',
   fr={'Arrhenius', A=2.89e13, n=0.0, T_a=-4.970e2*S},
   label='r13',
   note='identical to Stanford, r5'
}

reaction{
   'HO2 + HO2 <=> H2O2 + O2',
   fr={'Arrhenius', A=4.20e14, n=0.0, T_a=1.200e4*S},
   label='r14a',
   note='identical to Stanford, r6b'
}
reaction{
   'HO2 + HO2 <=> H2O2 + O2',
   fr={'Arrhenius', A=1.30e11, n=0.0, T_a=-1.630e3*S},
   label='r14b',
   note='identical to Stanford, r6a'
}

reaction{
   'H2O2 (+ M) <=> OH + OH (+ M)',
   fr={'pressure dependent',
       k_inf={A=2.00e12, n=0.90, T_a=4.875e4*S},
       k_0={A=2.49e24, n=-2.30, T_a=4.875e4*S},
       Troe={F_cent=0.42}
   },
   efficiencies={H2O=7.5,H2O2=7.7,O2=1.2,H2=3.7},
   label='r15',
   note='differs from Stanford, r3'
}

reaction{
   'H2O2 + H <=> H2O + OH',
   fr={'Arrhenius', A=2.41e13, n=0.0, T_a=3.970e3*S},
   label='r16',
   note='differs from Stanford, r16'
}

reaction{
   'H2O2 + H <=> HO2 + H2',
   fr={'Arrhenius', A=4.82e13, n=0.0, T_a=7.950e3*S},
   label='r17',
   note='differs from Stanford, r15'
}

reaction{
   'H2O2 + O <=> OH + HO2',
   fr={'Arrhenius', A=9.55e6, n=2.0, T_a=3.970e3*S},
   label='r18',
   note='differs from Stanford, r17'
}

reaction{
   'H2O2 + OH <=> HO2 + H2O',
   fr={'Arrhenius', A=1.74e12, n=0.0, T_a=3.180e2*S},
   label='r19a',
   note='identical to Stanford, r4a'
}
reaction{
   'H2O2 + OH <=> HO2 + H2O',
   fr={'Arrhenius', A=7.59e13, n=0.0, T_a=7.270e3*S},
   label='r19b',
   note='identical to Stanford, r4b'
}
