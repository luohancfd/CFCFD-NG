# Goecken_04.py
#
# This file provides a chemical kinetic description
# of Titan atmospheric entry,
#
# The reaction rates are taken from:
#
# Goecken (2004)
# N2-CH4-Ar Chemical Kinetic Model for Simulations of
# Atmospheric Entry to Titan
#
# This file prepared by...
# Rowan J Gollan
# 09-May-2006


elem_list = [ 'H', 'C', 'N' ]
ndata.declare_elements( elem_list )

species_list = [ 'N2', 'CH4', 'CH3', 'CH2', 'CH', 'C2', 'H2',
                 'CN', 'NH', 'HCN', 'N', 'C', 'H', 'N2+',
                 'CN+', 'N+', 'C+', 'H+', 'e-' ]
                 
                 
ndata.declare_species( species_list )

#
# Dissociation reactions
#

r1 = make_reaction("N2 + M <=> N + N + M",
                   fr=GA_forward( 7.00e21, -1.60, 113200.0 ),
                   efficiencies=[ ('N', 4.2857), ('C', 4.2857),
                                  ('H', 4.2857), ('e-', 428.57)  ] )

r2 = make_reaction("CH4 + M <=> CH3 + H + M",
                   fr=GA_forward( 4.70e47, -8.20, 59200.0 ) )

r3 = make_reaction("CH3 + M <=> CH2 + H + M",
                   fr=GA_forward( 1.02e16, 0.0, 45600.0 ) )

r4 = make_reaction("CH3 + M <=> CH + H2 + M",
                   fr=GA_forward( 5.00e15, 0.0, 42800.0 ) )

r5 = make_reaction("CH2 + M <=> CH + H + M",
                   fr=GA_forward( 4.00e15, 0.0, 41800.0 ) )

r6 = make_reaction("CH2 + M <=> C + H2 + M",
                   fr=GA_forward( 1.30e14, 0.0, 29700.0 ) )

r7 = make_reaction("CH + M <=> C + H + M",
                   fr=GA_forward( 1.90e14, 0.0, 33700.0 ) )

r8 = make_reaction("C2 + M <=> C + C + M",
                   fr=GA_forward( 1.50e16, 0.0, 716000.0 ) )

r9 = make_reaction("H2 + M <=> H + H + M",
                   fr=GA_forward( 2.23e14, 0.0, 48350.0 ) )

r10 = make_reaction("CN + M <=> C + N + M",
                    fr=GA_forward( 2.53e14, 0.0, 71000.0 ) )

r11 = make_reaction("NH + M <=> N + H + M",
                    fr=GA_forward( 1.80e14, 0.0, 37600.0 ) )

r12 = make_reaction("HCN + M <=> CN + H + M",
                    fr=GA_forward( 3.57e26, -2.60, 62845.0 ) )

#
# Radical reactions
# 

r13 = make_reaction("CH3 + N <=> HCN + H + H",
                    fr=GA_forward( 7.00e13, 0.00, 0.00 ) )

r14 = make_reaction("CH3 + H <=> CH2 + H2",
                    fr=GA_forward( 6.03e13, 0.00, 7600.0 ) )

r15 = make_reaction("CH2 + N2 <=> HCN + NH",
                    fr=GA_forward( 4.82e12, 0.0, 18000.0 ) )

r16 = make_reaction("CH2 + N <=> HCN + H",
                    fr=GA_forward( 5.00e13, 0.0, 0.0 ) )

r17 = make_reaction("CH2 + H <=> CH + H2",
                    fr=GA_forward( 6.03e12, 0.0, -900.0 ) )

r18 = make_reaction("CH + N2 <=> HCN + N",
                    fr=GA_forward( 4.40e12, 0.0, 11060.0 ) )

r19 = make_reaction("CH + C <=> C2 + H",
                    fr=GA_forward( 2.00e14, 0.00, 0.00 ) )

r20 = make_reaction("C2 + N2 <=> CN + CN",
                    fr=GA_forward( 1.50e13, 0.0, 21000.0 ) )

r21 = make_reaction("CN + H2 <=> HCN + H",
                    fr=GA_forward( 2.95e5, 0.0, 1130.0 ) )

r22 = make_reaction("CN + C <=> C2 + N",
                    fr=GA_forward( 5.00e13, 0.0, 13000.0 ) )

r23 = make_reaction("N + H2 <=> NH + H",
                    fr=GA_forward( 1.60e14, 0.0, 12650.0 ) )

r24 = make_reaction("C + N2 <=> CN + N",
                    fr=GA_forward( 5.24e13, 0.0, 22600.0 ) )

r25 = make_reaction("C + H2 <=> CH + H",
                    fr=GA_forward( 4.00e14, 0.0, 11700.0 ) )

r26 = make_reaction("H + N2 <=> NH + N",
                    fr=GA_forward( 3.00e12, 0.50, 71400.0 ) )

r27 = make_reaction("H + CH4 <=> CH3 + H2",
                    fr=GA_forward( 1.32e4, 3.00, 4045.0 ) )

#
# Ionization reactions
#

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

reactions_list = [ r1,  r2,  r3,  r4,  r5,
                   r6,  r7,  r8,  r9,  r10,
                   r11, r12, r13, r14, r15,
                   r16, r17, r18, r19, r20,
                   r21, r22, r23, r24, r25,
                   r26, r27, r28, r29, r30,
                   r31, r32, r34, r35 ]

ndata.declare_reactions( reactions_list )



                   




