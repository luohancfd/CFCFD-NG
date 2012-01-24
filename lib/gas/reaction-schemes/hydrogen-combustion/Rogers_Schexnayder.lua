-- Author: Rowan J. Gollan
-- Date: 29-Mar-2009
-- Place: Poquoson, Virginia, USA
--
-- Adapted from Python file: rogers_schexnayder.py
--
-- This file provides a reaction scheme for
-- hydrogen combustion in air.
-- NOTE: This scheme does not include carbonaceous compounds
-- or Argon (or any of the associated reactions).
--
-- Reference:
-- Rogers, R.C. and Schexnayder, Jr., C.J. (1981)
-- Chemical Kinetic Analysis of Hydroden-Air
-- Ignition and Reaction Times
-- NASA Technical Paper 1856
--
-- Species used: O, O2, N, N2, H, H2, H2O, HO2, OH, NO, NO2, HNO2, HNO3, O3, H2O2, HNO

reaction{
   'O2 + M <=> O + O + M',
   fr={"Arrhenius", A=0.72e19, n=-1.0, T_a=59340.0},
   efficiencies={O2=4.0, O=10.0, H2O=2.0},
   label='r1'
}

reaction{
   'M + H2 <=> H + H + M',
   fr={'Arrhenius', A=0.55e19, n=-1.0, T_a=51987.0},
   efficiencies={H=5.0, H2=2.0, H2O=8.0},
   label='r2'
}

reaction{
   'M + H2O <=> H + OH + M',
   fr={'Arrhenius', A=0.52e22, n=-1.5, T_a=59386.0},
   efficiencies={H2O=6.0},
   label='r3'
}

reaction{
   'H + O2 + M <=> HO2 + M',
   fr={'Arrhenius', A=0.23e16, n=0.0, T_a=-403.0},
   efficiencies={H2=2.0, H2O=13.0},
   label='r4'
}

reaction{
   'M + NO2 <=> NO + O + M',
   fr={'Arrhenius', A=0.11e17, n=0.0, T_a=32710.0},
   label='r5'
}

reaction{
   'M + NO <=> N + O + M',
   fr={'Arrhenius', A=0.41e19, n=-1.0, T_a=75330.0},
   label='r6'
}

-- not included, ignoring the carbonaceous compounds
r7 = 'M + O + CO <=> CO2 + M'

reaction{
   'M + H + NO <=> HNO + M',
   fr={'Arrhenius', A=0.54e16, n=0.0, T_a=-300.0},
   efficiencies={H2O=3.0},
   label='r8'
}

reaction{
   'M + H2O2 <=> OH + OH + M',
   fr={'Arrhenius', A=0.12e18, n=0.0, T_a=22899.0},
   efficiencies={H2O=6.0},
   label='r9'
}

reaction{
   'M + OH + NO <=> HNO2 + M',
   fr={'Arrhenius', A=0.80e16, n=0.0, T_a=-1000.0},
   label='r10'
}

reaction{
   'M + OH + NO2 <=> HNO3 + M',
   fr={'Arrhenius', A=0.13e17, n=0.0, T_a=-1107.0},
   label='r11'
}

reaction{
   'M + O3 <=> O2 + O + M',
   fr={'Arrhenius', A=0.13e22, n=-2.0, T_a=12800.0},
   efficiencies={O2=1.5},
   label='r12'
}

r13 = 'M + HCO <=> CO + H + M'

reaction{
   'M + O + H <=> OH + M',
   fr={'Arrhenius', A=0.71e19, n=-1.0, T_a=0.0},
   label='r14'
}

reaction{
   'H2O + O <=> OH + OH',
   fr={'Arrhenius', A=0.58e14, n=0.0, T_a=9059.0},
   label='r15'
}

reaction{
   'H2 + OH <=> H2O + H',
   fr={'Arrhenius', A=0.20e14, n=0.0, T_a=2600.0},
   label='r16'
}

reaction{
   'O2 + H <=> OH + O',
   fr={'Arrhenius', A=0.22e15, n=0.0, T_a=8455.0},
   label='r17'
}

reaction{
   'H2 + O <=> OH + H',
   fr={'Arrhenius', A=0.75e14, n=0.0, T_a=5586.0},
   label='r18'
}

reaction{
   'H2 + O2 <=> OH + OH',
   fr={'Arrhenius', A=0.10e14, n=0.0, T_a=21641.0},
   label='r19'
}

reaction{
   'H + HO2 <=> H2 + O2',
   fr={'Arrhenius', A=0.24e14, n=0.0, T_a=350.0},
   label='r20'
}

reaction{
   'H2 + O2 <=> H2O + O',
   fr={'Arrhenius', A=0.41e14, n=0.0, T_a=25400.0},
   label='r21'
}

reaction{
   'H + HO2 <=> OH + OH',
   fr={'Arrhenius', A=0.24e15, n=0.0, T_a=950.0},
   label='r22'
}

reaction{
   'H2O + O <=> H + HO2',
   fr={'Arrhenius', A=0.58e12, n=0.5, T_a=28686.0},
   label='r23'
}

reaction{
   'O + HO2 <=> OH + O2',
   fr={'Arrhenius', A=0.5e14, n=0.0, T_a=503.0},
   label='r24'
}

reaction{
   'OH + HO2 <=> O2 + H2O',
   fr={'Arrhenius', A=0.3e14, n=0.0, T_a=0.0},
   label='r25'
}

reaction{
   'H2 + HO2 <=> H2O + OH',
   fr={'Arrhenius', A=0.2e14, n=0.0, T_a=12582.0},
   label='r26'
}

reaction{
   'HO2 + H2 <=> H + H2O2',
   fr={'Arrhenius', A=0.73e12, n=0.0, T_a=9400.0},
   label='r27'
}

reaction{
   'H2O2 + H <=> OH + H2O',
   fr={'Arrhenius', A=0.32e15, n=0.0, T_a=4504.0},
   label='r28'
}

reaction{
   'HO2 + OH <=> O + H2O2',
   fr={'Arrhenius', A=0.52e11, n=0.5, T_a=10600.0},
   label='r29'
}

reaction{
   'HO2 + H2O <=> OH + H2O2',
   fr={'Arrhenius', A=0.28e14, n=0.0, T_a=16500.0},
   label='r30'
}

reaction{
   'HO2 + HO2 <=> H2O2 + O2',
   fr={'Arrhenius', A=0.2e13, n=0.0, T_a=0.0},
   label='r31'
}

reaction{
   'O + O3 <=> O2 + O2',
   fr={'Arrhenius', A=0.10e14, n=0.0, T_a=2411.0},
   label='r32'
}

reaction{
   'O3 + NO <=> NO2 + O2',
   fr={'Arrhenius', A=0.54e12, n=0.0, T_a=1200.0},
   label='r33'
}

reaction{
   'O3 + H <=> OH + O2',
   fr={'Arrhenius', A=0.70e14, n=0.0, T_a=560.0},
   label='r34'
}

reaction{
   'O3 + OH <=> O2 + HO2',
   fr={'Arrhenius', A=0.90e12, n=0.0, T_a=1000.0},
   label='r35'
}

reaction{
   'O + N2 <=> NO + N',
   fr={'Arrhenius', A=0.50e14, n=0.0, T_a=37940.0},
   label='r36'
}

reaction{
   'H + NO <=> OH + N',
   fr={'Arrhenius', A=0.17e15, n=0.0, T_a=24500.0},
   label='r37'
}

reaction{
   'O + NO <=> O2 + N',
   fr={'Arrhenius', A=0.15e10, n=1.0, T_a=19500.0},
   label='r38'
}

reaction{
   'NO2 + H <=> NO + OH',
   fr={'Arrhenius', A=0.35e15, n=0.0, T_a=740.0},
   label='r39'
}

reaction{
   'NO2 + O <=> NO + O2',
   fr={'Arrhenius', A=0.10e14, n=0.0, T_a=302.0},
   label='r40'
}

reaction{
   'NO2 + H2 <=> HNO2 + H',
   fr={'Arrhenius', A=0.24e14, n=0.0, T_a=14595.0},
   label='r41'
}

reaction{
   'HO2 + NO <=> NO2 + OH',
   fr={'Arrhenius', A=0.30e13, n=0.5, T_a=1208.0},
   label='r42'
}

reaction{
   'NO2 + H2O <=> HNO2 + OH',
   fr={'Arrhenius', A=0.32e13, n=0.0, T_a=22000.0},
   label='r43'
}

reaction{
   'NO2 + OH <=> HNO2 + O',
   fr={'Arrhenius', A=0.21e13, n=0.0, T_a=12580.0},
   label='r44'
}

r45 = 'CO + OH <=> CO2 + H'
r46 = 'CO2 + O <=> O2 + CO'
r47 = 'H2O + CO <=> HCO + OH'
r48 = 'OH + CO <=> HCO + O'
r49 = 'H2 + CO <=> HCO + H'
r50 = 'HO2 + CO <=> CO2 + OH'

reaction{
   'HNO + H <=> H2 + NO',
   fr={'Arrhenius', A=0.48e13, n=0.0, T_a=0.0},
   label='r51'
}

reaction{
   'HNO + OH <=> H2O + NO',
   fr={'Arrhenius', A=0.36e14, n=0.0, T_a=0.0},
   label='r52'
}

r53 = 'NO + CO <=> CO2 + N'
r54 = 'NO2 + CO <=> NO + CO2'

reaction{
   'NO + HO2 <=> HNO + O2',
   fr={'Arrhenius', A=0.72e13, n=0.5, T_a=5500.0},
   label='r55'
}

reaction{
   'HNO + O <=> NO + OH',
   fr={'Arrhenius', A=0.5e12, n=0.5, T_a=0.0},
   label='r56'
}

reaction{
   'HNO3 + O <=> HO2 + NO2',
   fr={'Arrhenius', A=0.10e12, n=0.0, T_a=0.0},
   label='r57'
}

reaction{
   'HO2 + NO2 <=> HNO2 + O2',
   fr={'Arrhenius', A=0.20e12, n=0.0, T_a=0.0},
   label='r58'
}

r59 = 'HCO + O2 <=> CO + HO2'

reaction{
   'O3 + HO2 <=> 2 O2 + OH',
   fr={'Arrhenius', A=0.10e12, n=0.0, T_a=1409.0},
   label='r60'
}


