-- single-step-reaction.lua
-- a simple fuel-air combustion mechanism.
-- PJ 22-Apr-2011

reaction{
   'AA + BB => 2AB',
   fr={'Arrhenius', A=1.0e6, n=0.0, T_a=2500.0},
}

