# nitrogen-5sp-6r-noneq.py
#
# The chemical kinetic system listed here provides
# a description of the simple dissociation
# of nitrogen.
#
# This file created by...
# Rowan J Gollan
# 
# And extended to include ionic reactions by...
# Daniel Potter
#
# 07-Mar-2006 : initial version
# 15-Mar-2006 : added K_eq options
# 11-Mar-2007 : added ionic reactions
# 01-Apr-2007 : amended for thermochemical nonequilibrium
#
# Sources: Park (1988) 'Two-temperature kinetic model for nitrogen'
#          Gokcen (2004) 'N2-CH4-Ar' Chemical kinetic model for simulations...'
#          CNES TC2T1 Reaction scheme
#
# Note:
#      -> Park's rates seem wrong, units might be not correct for some constants
#         eg 2.50e10 as opposed to 2.50e34 for A in r5
#      -> Gokcen's rates are a compilation of Park's work, and they almost match
#         CEA equilibrium conditions, so going with them for now.

ndata.method = "ode_mc"

USING_PARK = False
USING_GOKCEN = True
ndata.declare_elements( [ "e-", "N" ] )
ndata.declare_species( [  "N2", "N", "N2+", "N+", "e-" ] )
if USING_PARK:
    r0 = make_reaction("N2 + N2 <=> N + N + N2",
                       fr=GA_forward( 3.7e21, -1.6, 113200.0, 0, 1.0, 113200.0/3.0),
                       reac_type="Dissociation",
                       participating_ivib=[0])
    r1 = make_reaction( "N2 + N <=> N + N + N",
                        fr=ND_forward( 1.66e22, -1.6, 113200.0, 0, 1.0, 113200.0/3.0),
                        reac_type="Dissociation",
                        participating_ivib=[0])
    r2 = make_reaction( "N2 + e- <=> N + N + e-",
                        fr=ND_forward( 8.3e24, -1.6, 113200.0, 0, 1.0, 113200.0/3.0),
                        reac_type="Dissociation",
                        participating_ivib=[0])
    r3 = make_reaction( "N + N <=> N2+ + e-",
                        fr=GA_forward( 1.79e11, -0.77, 67500.0 ) )
    r4 = make_reaction( "N2 + N+ <=> N2+ + N",
    		        fr=GA_forward(9.85e12 , -0.18, 12100.0) )
    r5 = make_reaction("N + e- <=> N+ + e- + e-",
	        	fr=GA_forward(2.50e10 , -3.82, 168600.0) )

if USING_GOKCEN:
    r0 = make_reaction("N2 + N2 <=> N + N + N2",
                       fr=ND_forward( 7.0e21, -1.6, 113200.0, 0, 1.0, 113200.0/3.0),
                       reac_type="Dissociation",
                       participating_ivib=[0])
    r1 = make_reaction( "N2 + N <=> N + N + N",
                        fr=ND_forward( 3.0e22, -1.6, 113200.0, 0, 1.0, 113200.0/3.0),
                        reac_type="Dissociation",
                        participating_ivib=[0])
    r2 = make_reaction( "N2 + e- <=> N + N + e-",
                        fr=ND_forward( 3.0e24, -1.6, 113200.0, 0, 1.0, 113200.0/3.0),
                        reac_type="Dissociation",
                        participating_ivib=[0])
    r3 = make_reaction( "N + N <=> N2+ + e-",
                        fr=GA_forward( 4.40e7, 1.5, 67500.0 ) )
    r4 = make_reaction( "N2 + N+ <=> N2+ + N",
    		        fr=GA_forward(1.0e12 , 0.5, 12200.0) )
    r5 = make_reaction("N + e- <=> N+ + e- + e-",
	        	fr=GA_forward(2.50e34 , -3.82, 168600.0) )

ndata.declare_reactions( [ r0, r1, r2, r3, r4, r5 ] )
ndata.notes ="""#
# The full thermal dissociation and ionisation of molecular nitrogen.
# This set includes the effect of thermal nonequilibrium on the dissociation
# rate of nitrogen.
#"""
