-- Author: Rowan J. Gollan
-- Date: 11-Nov-2009
-- Place: Poquoson, Virginia, USA
--
-- Updated from the work by RJG and DFP as found in
-- lib/gas_models2/input_files/nitrogen/nitrogen-5sp-6r.py
--
-- Note: Based on Dan's comments, I've only included
-- the Goekcen rates at present.
--
-- Modified: Luke Doherty
-- Date: 09-Jun-2011
-- Place: UQ
-- Modifications:
--     : labelled each reaction
--     : introduced the WITH_IONIZATION flag following 
--     the work of F.Zander on the gupta_etal_air_reaction
--     scheme. This flag was introduced to allow for useage
--     with NENZFr.
--
-- Reference:
-- Goekcen (2004)
-- N2-CH4-Ar Chemical Kinetic Model for Simulations of 
-- Atmospheric Entry to Titan
-- AIAA Paper 2004-2469
--

WITH_IONIZATION = false

reaction{
   'N2 + N2 <=> N + N + N2',
   fr={'Arrhenius', A=7.0e21, n=-1.6, T_a=113200.0},
   label='r1'
}

reaction{
   'N2 + N <=> N + N + N',
   fr={'Arrhenius', A=3.0e22, n=-1.6, T_a=113200.0},
   label='r2'
}

reaction{
   'N2 + e- <=> N + N + e-',
   fr={'Arrhenius', A=3.0e24, n=-1.6, T_a=113200.0},
   label='r3'
}

reaction{
   'N + N <=> N2+ + e-',
   fr={'Arrhenius', A=4.40e7, n=1.5, T_a=67500.0},
   label='r4'
}

reaction{
   'N2 + N+ <=> N2+ + N',
   fr={'Arrhenius', A=1.0e12, n=0.5, T_a=12200.0},
   label='r5'
}

reaction{
   'N + e- <=> N+ + e- + e-',
   fr={'Arrhenius', A=2.50e34, n=-3.82, T_a=168600.0},
   label='r6'
}


if not WITH_IONIZATION then
	select_reactions_by_label({'r1', 'r2'})
end



