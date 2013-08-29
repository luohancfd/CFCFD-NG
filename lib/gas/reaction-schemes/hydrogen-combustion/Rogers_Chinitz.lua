-- Author: Wilson Chan
-- Date: 07-Oct-2010
-- Place: Centre for Hypersonics, UQ
--
-- Updates, tweaks, etc. by Fabian Zander and Rowan Gollan.
--
-- References
-- ----------
-- Rogers, R.C. and Chinitz, W. (1983)
-- Using a Global Hydrpgen-Air Combustion Model in Turbulent
-- Reacting Flow
-- AIAA Journal, 21:4, pp. 586--592
--
-- Berglund, M., Fedina, E., Fureby, C., and Tegner, J. (2010)
-- Finite Rate Chemistry Large-Eddy Simulation of Self-Ignition in a
-- Supersonic Combustion Ramjet
-- AIAA Journal, 48:3, March 2010.
-- 
-- Notes
-- -----
-- 1. The pre-exponential factor 'A' in the Rogers and Chinitz
-- rate is dependent on the equivalence ratio. Hence, the
-- global equivalence ratio should be set when using this model.
--
-- 2. I (Rowan Gollan) don't recommend this scheme. It's troublesome
-- in a numerical sense. Rogers and Chinitz discuss some of this
-- in their paper.
--

phi = 1.0
A4 = (8.917*phi + 31.433/phi - 28.950)*1e47
A5 = (2.000 + 1.333/phi - 0.833*phi)*1e64

R_cal = 1.9872041 -- cal/(K.mol)

reaction{
   'H2 + O2 <=> 2OH',
   fr={'Arrhenius', A=A4, n=-10, T_a=4865/R_cal},
   label='r1'
}

reaction{
   '2OH + H2 <=> 2H2O',
   fr={'Arrhenius', A=A5, n=-13, T_a=42500/R_cal},
   label='r2'
}


