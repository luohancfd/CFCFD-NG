-- Author: Rowan J. Gollan
-- Date: 2015-03-17
--
-- Reference:
-- Hong, Z., Davidson, D.F. and Hanson, R.K. (2011)
-- An improved H2/O2 mechanism based on recent
-- shock tube/laser absorption measurements.
-- Combustion and Flame, 158, pp. 633--644
--
-- NOTE:
-- Table 1 in Hong et al contains the suggested
-- reaction mechanism with reaction rates.
-- However, there is also a Chemkin input file
-- available online where the article is hosted.
-- The input file has some differences with regards
-- to efficiency values. I have adopted those
-- values from the supplied input file.

S = 1.0/1.987

reaction{
   'H + O2 <=> OH + O',
   fr={'Arrhenius', A=1.04e14, n=0.0, T_a=15286*S},
   label='r1'
}

reaction{
   'H + O2 (+ M) <=> HO2 (+ M)',
   fr={'pressure dependent',
        k_inf={A=5.59e13, n=0.2, T_a=0.0},
	k_0={A=3.70e19, n=-1.0, T_a=0.0},
	Troe={F_cent=0.8}
   },
   efficiencies={H2O=1.0,H=0.0,O2=0.0,OH=0.0,O=0.0,HO2=0.0,H2O2=0.0},
   label='r2b'
}

reaction{
   'H + O2 (+ M) <=> HO2 (+ M)',
   fr={'pressure dependent',
        k_inf={A=5.59e13, n=0.2, T_a=0.0},
        k_0={A=5.69e18, n=-1.1, T_a=0.0},
        Troe={F_cent=0.7}
   },
   efficiencies={O2=1.0,H2O=0.0,H=0.0,OH=0.0,O=0.0,HO2=0.0,H2O2=0.0},
   label='r2c'
}
     
reaction{
   'H + O2 (+ M) <=> HO2 (+ M)',
   fr={'pressure dependent',
        k_inf={A=5.59e13, n=0.2, T_a=0.0},
	k_0={A=2.65e19, n=-1.3, T_a=0.0},
        Troe={F_cent=0.7}
   },
   efficiencies={H2=2.5,H2O2=12.0,H2O=0.0,O2=0.0},
   label='r2d'
}

reaction{
   'H2O2 (+ M) <=> 2OH (+ M)',
   fr={'pressure dependent',
        k_inf={A=8.59e14, n=0.0, T_a=48560*S},
        k_0={A=9.55e15, n=0.0, T_a=42203*S},
        Troe={F_cent=1.0}
   },
   efficiencies={N2=1.5,H2=2.5,H2O=15,H2O2=15},
   label='r3'
}

-- Reaction 4 appears twice.
-- The reaction rate constants are added together
-- as proposed by Hong et al in Section 2.5.
-- To achieve the same effect, the reaction can just
-- be listed twice with different reaction rates.
reaction{
   'OH + H2O2 <=> H2O + HO2',
   fr={'Arrhenius', A=1.74e12, n=0.0, T_a=318*S},
   label='r4a'
}
reaction{
   'OH + H2O2 <=> H2O + HO2',
   fr={'Arrhenius', A=7.59e13, n=0.0, T_a=7269*S},
   label='r4b'
}

reaction{
   'OH + HO2 <=> H2O + O2',
   fr={'Arrhenius', A=2.89e13, n=0.0, T_a=-500*S},
   label='r5'
}

reaction{
   'HO2 + HO2 <=> H2O2 + O2',
   fr={'Arrhenius', A=1.30e11, n=0.0, T_a=-1603*S},
   label='r6a'
}
reaction{
   'HO2 + HO2 <=> H2O2 + O2',
   fr={'Arrhenius', A=4.20e14, n=0.0, T_a=11980*S},
   label='r6b'
}

reaction{
   'H2O + M <=> H + OH + M',
   fr={'Arrhenius', A=6.06e27, n=-3.31, T_a=120770*S},
   efficiencies={O2=1.5,H2=3.0,H2O=0.0},
   label='r7a'
}
reaction{
   'H2O + H2O <=> H + OH + H2O',
   fr={'Arrhenius', A=1.00e26, n=-2.44, T_a=120160*S},
   label='r7b'
}

reaction{
   'OH + OH <=> H2O + O',
   fr={'Arrhenius', A=3.57e4, n=2.4, T_a=-2111*S},
   label='r8'
}

reaction{
   'O + H2 <=> H + OH',
   fr={'Arrhenius', A=3.82e12, n=0.0, T_a=7948*S},
   label='r9a'
}
reaction{
   'O + H2 <=> H + OH',
   fr={'Arrhenius', A=8.79e14, n=0.0, T_a=19170*S},
   label='r9b'
}

reaction{
   'H2 + OH <=> H2O + H',
   fr={'Arrhenius', A=2.17e8, n=1.52, T_a=3457*S},
   label='r10'
}

reaction{
   'H + HO2 <=> OH + OH',
   fr={'Arrhenius', A=7.08e13, n=0.0, T_a=300*S},
   label='r11'
}

reaction{
   'H + HO2 <=> H2O + O',
   fr={'Arrhenius', A=1.45e12, n=0.0, T_a=0.0},
   label='r12'
}

reaction{
   'H + HO2 <=> H2 + O2',
   fr={'Arrhenius', A=3.66e6, n=2.087, T_a=-1450*S},
   label='r13'
}

reaction{
   'O + HO2 <=> OH + O2',
   fr={'Arrhenius', A=1.63e13, n=0.0, T_a=-445*S},
   label='r14'
}

reaction{
   'H2O2 + H <=> HO2 + H2',
   fr={'Arrhenius', A=1.21e7, n=2.0, T_a=5200*S},
   label='r15'
}

reaction{
   'H2O2 + H <=> H2O + OH',
   fr={'Arrhenius', A=1.02e13, n=0.0, T_a=3577*S},
   label='r16'
}

reaction{
   'H2O2 + O <=> OH + HO2',
   fr={'Arrhenius', A=8.43e11, n=0.0, T_a=3970*S},
   label='r17'
}

reaction{
   'H2 + M <=> H + H + M',
   fr={'Arrhenius', A=5.84e18, n=-1.1, T_a=104380*S},
   efficiencies={H20=14.4,H2O2=14.4,H2=0.0,O2=0.0},
   label='r18a'
}
reaction{
   'H2 + H2 <=> H + H + H2',
   fr={'Arrhenius', A=9.03e14, n=0.0, T_a=96070*S},
   label='r18b'
}
reaction{
   'H2 + O2 <=> H + H + O2',
   fr={'Arrhenius', A=4.58e19, n=-1.4, T_a=104380*S},
   label='r18c'
}

reaction{
   'O + O + M <=> O2 + M',
   fr={'Arrhenius', A=6.16e15, n=-0.5, T_a=0.0},
   efficiencies={H2=2.5,H2O=12,H2O2=12},
   label='r19'
}

reaction{
   'O + H + M <=> OH + M',
   fr={'Arrhenius', A=4.71e18, n=-1.0, T_a=0.0},
   efficiencies={H2=2.5,H2O=12,H2O2=12},
   label='r20'
}




