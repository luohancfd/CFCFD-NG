-- Author: Rowan J. Gollan
-- Date: 09-Jun-2009
-- Place: Poquoson, Virginia, USA
--
-- This files file provides a chemical kinetic description
-- of Titan atmospheric entry,
--
-- The reaction rates are taken from:
--
-- Goekcen (2004)
-- N2-CH4-Ar Chemical Kinetic Model for Simulations of
-- Atmospheric Entry to Titan
--
-- This is adapted from the the Python file: Goekcen_04.py
--
--

reaction{
   'N2 + M <=> N + N + M',
   fr={'Arrhenius', A=7.00e21, n=-1.60, T_a=113200.0},
--   efficiencies={N=4.2857, C=4.2857, H=4.2857, e-=428.57}, -- DON'T include 'e-' presently
   efficiencies={N=4.2857, C=4.2857, H=4.2857},
   label='r1'
}

reaction{
   'CH4 + M <=> CH3 + H + M',
   fr={'Arrhenius', A=4.70e47, n=-8.20, T_a=59200.0},
   label='r2'
}

reaction{
   'CH3 + M <=> CH2 + H + M',
   fr={'Arrhenius', A=1.02e16, n=0.0, T_a=45600.0},
   label='r3'
}

reaction{
   'CH3 + M <=> CH + H2 + M',
   fr={'Arrhenius', A=5.00e15, n=0.0, T_a=42800.0},
   label='r4'
}

reaction{
   'CH2 + M <=> CH + H + M',
   fr={'Arrhenius', A=4.00e15, n=0.0, T_a=41800.0},
   label='r5'
}

reaction{
   'CH2 + M <=> C + H2 + M',
   fr={'Arrhenius', A=1.30e14, n=0.0, T_a=29700.0},
   label='r6'
}

reaction{
   'CH + M <=> C + H + M',
   fr={'Arrhenius', A=1.90e14, n=0.0, T_a=33700.0},
   label='r7'
}

reaction{
   'C2 + M <=> C + C + M',
   fr={'Arrhenius', A=1.50e16, n=0.0, T_a=716000.0},
   label='r8'
}

reaction{
   'H2 + M <=> H + H + M',
   fr={'Arrhenius', A=2.23e14, n=0.0, T_a=48350.0},
   label='r9'
}

reaction{
   'CN + M <=> C + N + M',
   fr={'Arrhenius', A=2.53e14, n=0.0, T_a=71000.0},
   label='r10'
}

reaction{
   'NH + M <=> N + H + M',
   fr={'Arrhenius', A=1.80e14, n=0.0, T_a=37600.0},
   label='r11'
}

reaction{
   'HCN + M <=> CN + H + M',
   fr={'Arrhenius', A=3.57e26, n=-2.60, T_a=62845.0},
   label='r12'
}

------           ------
-- Radical reactions --
------           ------ 

reaction{
   'CH3 + N <=> HCN + H + H',
   fr={'Arrhenius', A=7.00e13, n=0.00, T_a=0.00},
   label='r13'
}

reaction{
   'CH3 + H <=> CH2 + H2',
   fr={'Arrhenius', A=6.03e13, n=0.00, T_a=7600.0},
   label='r14'
}

reaction{
   'CH2 + N2 <=> HCN + NH',
   fr={'Arrhenius', A=4.82e12, n=0.0, T_a=18000.0},
   label='r15'
}

reaction{
   'CH2 + N <=> HCN + H',
   fr={'Arrhenius', A=5.00e13, n=0.0, T_a=0.0},
   label='r16'
}

reaction{
   'CH2 + H <=> CH + H2',
   fr={'Arrhenius', A=6.03e12, n=0.0, T_a=-900.0},
   label='r17'
}

reaction{
   'CH + N2 <=> HCN + N',
   fr={'Arrhenius', A=4.40e12, n=0.0, T_a=11060.0},
   label='r18'
}

reaction{
   'CH + C <=> C2 + H',
   fr={'Arrhenius', A=2.00e14, n=0.00, T_a=0.00},
   label='r19'
}

reaction{
   'C2 + N2 <=> CN + CN',
   fr={'Arrhenius', A=1.50e13, n=0.0, T_a=21000.0},
   label='r20'
}

reaction{
   'CN + H2 <=> HCN + H',
   fr={'Arrhenius', A=2.95e5, n=0.0, T_a=1130.0},
   label='r21'
}

reaction{
   'CN + C <=> C2 + N',
   fr={'Arrhenius', A=5.00e13, n=0.0, T_a=13000.0},
   label='r22'
}

reaction{
   'N + H2 <=> NH + H',
   fr={'Arrhenius', A=1.60e14, n=0.0, T_a=12650.0},
   label='r23'
}

reaction{
   'C + N2 <=> CN + N',
   fr={'Arrhenius', A=5.24e13, n=0.0, T_a=22600.0},
   label='r24'
}

reaction{
   'C + H2 <=> CH + H',
   fr={'Arrhenius', A=4.00e14, n=0.0, T_a=11700.0},
   label='r25'
}

reaction{
   'H + N2 <=> NH + N',
   fr={'Arrhenius', A=3.00e12, n=0.50, T_a=71400.0},
   label='r26'
}

reaction{
   'H + CH4 <=> CH3 + H2',
   fr={'Arrhenius', A=1.32e4, n=3.00, T_a=4045.0},
   label='r27'
}

--[=[ NOT YET INCLUDED.
#
# Ionization reactions
#

if WITH_IONIZATION:

    r28 = make_reaction("N + N <=> N2+ + e-",
                        fr=GA_forward( 4.40e7, 1.50, 67500.0 ) )

    r29 = make_reaction("C + N <=> CN+ + e-",
                        fr=GA_forward( 1.00e15, 1.50, 164400.0 ) )

    r30 = make_reaction("N + e- <=> N+ + e- + e-",
                        fr=GA_forward( 2.50e34, -3.82, 168600.0 ) )

    r31 = make_reaction("C + e- <=> C+ + e- + e-",
                        fr=GA_forward( 3.70e31, -3.00, 130720.0 ) )

    r32 = make_reaction("H + e- <=> H+ + e- + e-",
                        fr=GA_forward( 2.20e30, -2.80, 157800.0 ) )

    # r33 omitted : Argon reaction

    r34 = make_reaction("CN+ + N <=> CN + N+",
                        fr=GA_forward( 9.80e12, 0.00, 40700.0 ) )

    r35 = make_reaction("C+ + N2 <=> N2+ + C",
                        fr=GA_forward( 1.11e14, -0.11, 50000.0 ) )

--]=]