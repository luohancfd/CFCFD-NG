-- Author: Wilson Chan
-- Date: 07-Oct-2010
-- Place: Centre for Hypersonics, UQ
--
-- This file provides one chemical kinetic descriptions
-- of hydrogen combustion.  You can select between the various
-- options below by setting the 'model' variable below to one of
-- the strings listed below.
--
-- IN_AIR   : a 4-species, 2-reactions description of hydrogen
--            combustion in air (N2 and O2)
--
-- The numbering of reactions in this file corresponds to
-- Table 2 in Berglund et al (2010).
--
-- Reference:
--  Evans, J.S. and Shexnayder Jr, C.J. (1980)
--
-- Berglund, M., Fedina, E., Fureby, C., and Tegner, J.
-- ``Finite Rate Chemistry Large-Eddy Simulation of Self-Ignition in a
-- Supersonic Combustion Ramjet,'' AIAA Journal, Vol. 48, No. 3, March 2010.
--
-- Updated: FZ, 18th Oct, St Lucia
-- Arrhenius parameters changed to reflect the correct values in Berglund paper
--
-- Species used: H2, O2, OH, H2O

reaction{
   'H2 + O2 <=> OH + OH',
   fr={'Arrhenius', A=2.30e16, n=0.0, T_a=5134.0},
   label='r1'
}

reaction{
   'OH + OH + H2 <=> H2O + H2O',
   fr={'Arrhenius', A=1.83e18, n=0.0, T_a=11067.0},
   label='r2'
}