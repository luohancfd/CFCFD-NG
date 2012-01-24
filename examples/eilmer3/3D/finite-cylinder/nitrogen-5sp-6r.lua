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
-- Reference:
-- Goekcen (2004)
-- N2-CH4-Ar Chemical Kinetic Model for Simulations of 
-- Atmospheric Entry to Titan
-- AIAA Paper 2004-2469
--

reaction{
   'N2 + N2 <=> N + N + N2',
   fr={'Arrhenius', A=7.0e21, n=-1.6, T_a=113200.0}
}

reaction{
   'N2 + N <=> N + N + N',
   fr={'Arrhenius', A=3.0e22, n=-1.6, T_a=113200.0}
}

reaction{
   'N2 + e- <=> N + N + e-',
   fr={'Arrhenius', A=3.0e24, n=-1.6, T_a=113200.0}
}

reaction{
   'N + N <=> N2+ + e-',
   fr={'Arrhenius', A=4.40e7, n=1.5, T_a=67500.0}
}

reaction{
   'N2 + N+ <=> N2+ + N',
   fr={'Arrhenius', A=1.0e12, n=0.5, T_a=12200.0}
}

reaction{
   'N + e- <=> N+ + e- + e-',
   fr={'Arrhenius', A=2.50e34, n=-3.82, T_a=168600.0}
}

