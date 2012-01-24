-- Author: Rowan J. Gollan
-- Date: 02-Feb-2010
-- Place: Poquoson, Virginia, USA
--
-- Adapted from Python file: evans_scheznayder.py
--
-- This file provides four chemical kinetic descriptions
-- of hydrogen combustion.  You can select between the various
-- options below by setting the 'model' variable below to one of
-- the strings listed below.
--
-- REDUCED  : a 7-species, 8-reactions description of hydrogen
--            combustion in pure oxygen
-- PURE_O2  : a 7-species, 16-reactions description of hydrogen
--            combustion in pure oxygen
-- IN_AIR   : a 12-species, 25-reactions description of hydrogen
--            combustion in air (N2 and O2)
-- INERT_N2 : an 8-species, 16-reactions description of hydrogen
--            combustion in air with inert N2 (acting as diluent only).
--
-- The numbering of reactions in this file corresponds to
-- Table 1 in Evans and Schexnayder (1980).
--
-- Reference:
--  Evans, J.S. and Shexnayder Jr, C.J. (1980)
--  Influence of Chemical Kinetics and Unmixedness
--  on Burning in Supersonic Hydrogen Flames
--  AIAA Journal 18:2 pp 188--193
--
-- History:
--  07-Mar-2006 -- first prepared
--

options = {
   REDUCED=true,
   PURE_O2=true,
   IN_AIR=true,
   INERT_N2=true
}

-- User selects model here
model = 'INERT_N2'

-- Check that selection is valid
if options[model] == nil then
   print("User selected model: ", model)
   print("is not valid.")
   print("Valid models are:")
   for m,_ in pairs(options) do
      print(m)
   end
end

reaction{
   'HNO2 + M <=> NO + OH + M',
   fr={'Arrhenius', A=5.0e17, n=-1.0, T_a=25000.0},
   br={'Arrhenius', A=8.0e15, n=0.0, T_a=-1000.0},
   label='r1'
}

reaction{
   'NO2 + M <=> NO + O + M',
   fr={'Arrhenius', A=1.1e16, n=0.0, T_a=32712.0},
   br={'Arrhenius', A=1.1e15, n=0.0, T_a=-941.0},
   label='r2'
}

reaction{
   'H2 + M <=> H + H + M',
   fr={'Arrhenius', A=5.5e18, n=-1.0, T_a=51987.0},
   br={'Arrhenius', A=1.8e18, n=-1.0, T_a=0.0},
   label='r3'
}

reaction{
   'O2 + M <=> O + O + M',
   fr={'Arrhenius', A=7.2e18, n=-1.0, T_a=59340.0},
   br={'Arrhenius', A=4.0e17, n=-1.0, T_a=0.0},
   label='r4'
}

reaction{
   'H2O + M <=> OH + H + M',
   fr={'Arrhenius', A=5.2e21, n=-1.5, T_a=59386.0},
   br={'Arrhenius', A=4.4e20, n=-1.5, T_a=0.0},
   label='r5'
}

reaction{
   'OH + M <=> O + H + M',
   fr={'Arrhenius', A=8.5e18, n=-1.0, T_a=50830.0},
   br={'Arrhenius', A=7.1e18, n=-1.0, T_a=0.0},
   label='r6'
}

reaction{
   'HO2 + M <=> H + O2 + M',
   fr={'Arrhenius', A=1.7e16, n=0.0, T_a=23100.0},
   br={'Arrhenius', A=1.1e16, n=0.0, T_a=-440.0},
   label='r7'
}

reaction{
   'H2O + O <=> OH + OH',
   fr={'Arrhenius', A=5.8e13, n=0.0, T_a=9059.0},
   br={'Arrhenius', A=5.3e12, n=0.0, T_a=503.0},
   label='r8'
}

reaction{
   'H2O + H <=> OH + H2',
   fr={'Arrhenius', A=8.4e13, n=0.0, T_a=10116.0},
   br={'Arrhenius', A=2.0e13, n=0.0, T_a=2600.0},
   label='r9'
}

reaction{
   'O2 + H <=> OH + O',
   fr={'Arrhenius', A=2.2e14, n=0.0, T_a=8455.0},
   br={'Arrhenius', A=1.5e13, n=0.0, T_a=0.0},
   label='r10'
}

reaction{
   'H2 + O <=> OH + H',
   fr={'Arrhenius', A=7.5e13, n=0.0, T_a=5586.0},
   br={'Arrhenius', A=3.0e13, n=0.0, T_a=4429.0},
   label='r11'
}

reaction{
   'H2 + O2 <=> OH + OH',
   fr={'Arrhenius', A=1.7e13, n=0.0, T_a=24232.0},
   br={'Arrhenius', A=5.7e11, n=0.0, T_a=14922.0},
   label='r12'
}

reaction{
   'H2 + O2 <=> H + HO2',
   fr={'Arrhenius', A=1.9e13, n=0.0, T_a=24100.0},
   br={'Arrhenius', A=1.3e13, n=0.0, T_a=0.0},
   label='r13'
}

reaction{
   'OH + OH <=> H + HO2',
   fr={'Arrhenius', A=1.7e11, n=0.5, T_a=21137.0},
   br={'Arrhenius', A=6.0e13, n=0.0, T_a=0.0},
   label='r14'
}

reaction{
   'H2O + O <=> H + HO2',
   fr={'Arrhenius', A=5.8e11, n=0.5, T_a=28686.0},
   br={'Arrhenius', A=3.0e13, n=0.0, T_a=0.0},
   label='r15'
}

reaction{
   'OH + O2 <=> O + HO2',
   fr={'Arrhenius', A=3.7e11, n=0.64, T_a=27840.0},
   br={'Arrhenius', A=1.0e13, n=0.0, T_a=0.0},
   label='r16'
}

reaction{
   'H2O + O2 <=> OH + HO2',
   fr={'Arrhenius', A=2.0e11, n=0.5, T_a=36296.0},
   br={'Arrhenius', A=1.2e13, n=0.0, T_a=0.0},
   label='r17'
}

reaction{
   'H2O + OH <=> H2 + HO2',
   fr={'Arrhenius', A=1.2e12, n=0.21, T_a=39815.0},
   br={'Arrhenius', A=1.7e13, n=0.0, T_a=12582.0},
   label='r18'
}

reaction{
   'O + N2 <=> N + NO',
   fr={'Arrhenius', A=5.0e13, n=0.0, T_a=37940.0},
   br={'Arrhenius', A=1.1e13, n=0.0, T_a=0.0},
   label='r19'
}

reaction{
   'H + NO <=> N + OH',
   fr={'Arrhenius', A=1.7e14, n=0.0, T_a=24500.0},
   br={'Arrhenius', A=4.5e13, n=0.0, T_a=0.0},
   label='r20'
}

reaction{
   'O + NO <=> N + O2',
   fr={'Arrhenius', A=2.4e11, n=0.5, T_a=19200.0},
   br={'Arrhenius', A=1.0e12, n=0.5, T_a=3120.0},
   label='r21'
}

reaction{
   'NO + OH <=> H + NO2',
   fr={'Arrhenius', A=2.0e11, n=0.5, T_a=15500.0},
   br={'Arrhenius', A=3.5e14, n=0.0, T_a=740.0},
   label='r22'
}


reaction{
   'NO + O2 <=> O + NO2',
   fr={'Arrhenius', A=1.0e12, n=0.0, T_a=22800.0},
   br={'Arrhenius', A=1.0e13, n=0.0, T_a=302.0},
   label='r23'
}

reaction{
   'NO2 + H2 <=> H + HNO2',
   fr={'Arrhenius', A=2.4e13, n=0.0, T_a=14500.0},
   br={'Arrhenius', A=5.0e11, n=0.5, T_a=1500.0},
   label='r24'
}

reaction{
   'NO2 + OH <=> NO + HO2',
   fr={'Arrhenius', A=1.0e11, n=0.5, T_a=6000.0},
   br={'Arrhenius', A=3.0e12, n=0.5, T_a=1200.0},
   label='r25'
}

reactions_list = {}

if model == 'REDUCED' then
   reactions_list = {'r3', 'r4', 'r5', 'r6', 'r8', 'r9', 'r10', 'r11'}
end

if model == 'PURE_O2' or model == 'INERT_N2' then
   reactions_list = {'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10',
		     'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18'}
end


if model ~= 'IN_AIR' then
   -- For all other models we select only a subset.
   select_reactions_by_label(reactions_list)
end
