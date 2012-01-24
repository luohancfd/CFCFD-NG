# air-7sp-10r.py
#
# This file provides a chemical kinetic description
# of air chemistry typical of re-entry type problems.
#
# The reaction rates are taken from Table II of..
#
# Marrone, P.V. (1963)
# Inviscid, nonequilibrium flow behind bow and normal shock
# waves, Part I. General analysis and numerical examples
# CAL Report No. QM-16260A012 (I)
#
# Note: the final reaction in Table II is excluded as
# Marrone did not consider this in his example for air at M=12.3
#
# This file prepared by...
# Rowan J Gollan
# 03-May-2006

MARRONE_SCHEME = False
GUPTA_SCHEME = True


elem_list = [ 'O', 'N' ]
ndata.declare_elements( elem_list )

species_list = [ 'O2', 'N2', 'NO', 'O', 'N', 'NO+', 'e-' ]
ndata.declare_species( species_list )


r1 = make_reaction("O2 + O <=> O + O + O",
                   fr=GA_forward( 2.1e18, -0.5, 59380.0 ),
                   reac_type="Dissociation")

r2 = make_reaction("O2 + O2 <=> O + O + O2",
                   fr=GA_forward( 3.6e21, -1.5, 59380.0 ),
                   reac_type="Dissociation")

r3 = make_reaction("O2 + M <=> O + O + M",
                   fr=GA_forward( 1.19e21, -1.5, 59380.0 ),
                   efficiencies=[ ("O", 0.0), ("O2", 0.0) ] )

r4 = make_reaction("N2 + N <=> N + N + N",
                   fr=kf_from_eq_const(),
                   br=GA_backward( 7.5e20, -1.5, 0.0 ),
                   reac_type="Dissociation" )

r5 = make_reaction("N2 + N2 <=> N + N + N2",
                   fr=kf_from_eq_const(),
                   br=GA_backward( 1.5e20, -1.5, 0.0 ),
                   reac_type="Dissociation" )

r6 = make_reaction("N2 + M <=> N + N + M",
                   fr=kf_from_eq_const(),
                   br=GA_backward( 5e19, -1.5, 0.0 ),
                   efficiencies=[ ("N", 0.0), ("N2", 0.0) ] )

r7 = make_reaction("NO + M <=> N + O  + M",
                   fr=GA_forward( 5.18e21, -1.5, 75490.0 ) )
    
r8 = make_reaction("N + O2 <=> NO + O",
                   fr=GA_forward( 1.0e12, 0.5, 3120.0 ) )

r9 = make_reaction("O + N2 <=> NO + N",
                   fr=GA_forward( 5.0e13, 0.0, 38016.0 ) )

r10 = make_reaction("N + O <=> NO+ + e-",
                    fr=kf_from_eq_const(),
                    br=GA_backward( 1.80e21, -1.5, 0.0, nc=2 ) )

if GUPTA_SCHEME:
    r1.set_fr_model( GA_forward(3.61e18, -1.0, 59400.0) )

    r2.set_fr_model( GA_forward(3.61e18, -1.0, 59400.0) )

    r3 = make_reaction("O2 + M <=> O + O + M",
                       fr=GA_forward(3.61e18, -1.0, 59400.0),
                       efficiencies=[ ("O", 0.0), ("O2", 0.0) ] )

    r4.set_fr_model( GA_forward(4.15e22, -1.5, 113100.0) )
    r4.set_br_model( kb_from_eq_const() )

    r5.set_fr_model( GA_forward(1.92e17, -0.5, 113100.0) )
    r5.set_br_model( kb_from_eq_const() )

    r6 = make_reaction("N2 + M <=> N + N + M",
                       fr=GA_forward(1.92e17, -0.5, 113100.0),
                       efficiencies=[ ("N", 0.0), ("N2", 0.0) ] )

    r7.set_fr_model( GA_forward(3.97e20, -1.5, 75600.0) )

    r8.set_fr_model( GA_forward(9.63e11, 0.5, 3600.0) )

    r9.set_fr_model( GA_forward(6.75e13, 0.0, 37500.0) )

    r10.set_fr_model( GA_forward(9.03e9, 0.5, 32400.0) )


reactions_list = [ r1, r2, r3, r4, r5, r6, r7, r8, r9, r10 ]
ndata.declare_reactions( reactions_list )

ndata.notes = """#
# This file provides a chemical kinetic description
# of air chemistry typical of re-entry type problems
#
# Reference:
# Marrone, P.V. (1963)
# Inviscid, nonequilibrium flow behind bow and normal shock
# waves, Part I. General analysis and numerical examples
# CAL Report No. QM-16260A012 (I)
#"""


                   




