# ParkX-NE.py
#
# CO2-N2 reaction scheme with base reactions from Park (1994), plus additional 
# and modified reactions from Gokcen (2004), Tsang (1991) and Losev (1999), as 
# described in Lee (2007).  Two N2+ forming reactions and two reactions involving
# the CN+ cation from Gokcen (2004) that are not considered by Lee (2007) are 
# also included.  This chemistry scheme is applied in a two temperature
# context (thermal nonequilibrium). VC coupling is via the rate controlling
# temperatures recommended by Park.  CV coupling is non-preferential (nominal)
# for the moment, with a slightly better method soon to come. References:
#
# C.S. Park, J.T. Howe, R.L. Jaffe, and G.V. Candler.
# Review of chemical-kinetic problems of future NASA missions, II: Mars entries.
# Journal of Thermophysics and Heat Transfer, 8(1):9 - 22, 1994.
#
# Tahir Gokcen.
# N2 - CH4 - Ar Chemical Kinetic Model for Simulations of Atmospheric Entry to Titan.
# 37th AIAA Thermophysics Conference, 28 June - 1 July 2004, Portland, Oregon.
#
# W. Tsang and J.T. Herron. 
# Chemical kinetic data base for propellant combustion I. Reactions involving 
# NO, NO2 , HNO, HNO2 , HCN and N2O.
# Journal of Physical and Chemical Reference Data, 20(4):609 - 663, 1991. 
#
# S.A. Losev and G.G. Chernyi.
# Development of thermal protection system for interplanetary flight.
# Report 036-96, International Science and Technology Center, 
# Moscow State University, August 1999.
#
# Eun-Seok Lee, Chul Park, and Keun-Shik Chang.
# Shock-tube determination of CN formation rate in a CO-N2 mixture.
# Journal of Thermophysics and Heat Transfer, 21(1):50 - 56, 2007.
#
# NOTE: Reactions r0 -> r26 are taken verbatim from Lee (2007) p55.
#       Charged reactions, r27 -> r39
#
# Log:
# 08-Mar-08: 
# Foundation
#
# 08-Apr-08
# Working combination (TTTT)
# Unstable Gokcen CN reaction 30 is omitted
#
# 26-May-08
# Addition of optional Rond and Gorelov reaction rates
#
# 23-Aug-08
# Addition of Treanor-Marrone coupling
#
# 17-Sep-08
# Elec coupling

ndata.T_trigger = 300.0
ndata.method = "ode_mc"
ndata.elec_coupling = True

elem_list = [ "Ar", "C", "N", "O", "e-" ]
ndata.declare_elements( elem_list )

species_list = [ "CO2",  "CO",  "N2",  "CN",  "NO", "O2", "C2", "NCO", 
                 "CO+", "N2+", "CN+", "NO+", "O2+", "Ar",  "C",   "N",
		   "O",  "C+",  "N+", "O+",  "e-" ]
		   
# parkX_tTg.inp:
# species  =  CO2  CO  N2  CN  NO  O2  C2  NCO  CO+  N2+ CN+  NO+  O2+  Ar  C   N   O   C+  N+  O+  e-
# index    =  0    1   2   3   4   5   6   7    8    9   10   11   12   13  14  15  16  17  18  19  20

ndata.declare_species( species_list )

# 1.0 Dissociation Reactions

# 	1.1 C2 Dissociation reaction

r1 = make_reaction("C2 + M <=> C + C + M",
    		fr=Park_NE_forward( 1.5e16, 0.0, 69900.0, 6, 0.5, 0.0 ),
                reac_type="ThirdBodyReaction",
		participating_ivib=[6] )

# 	1.2 N2 Dissociation reactions

eA = 3.000e22 / 7.000e21
eE = 3.000e24 / 7.000e21

r2 = make_reaction("N2 + M <=> N + N + M",
		fr=Park_NE_forward( 7.000e21, -1.6, 113200.0, 2, 0.5, 0.0 ),
                efficiencies=[ ("CO2", 1.0), ( "CO", 1.0), ( "N2", 1.0), ( "CN", 1.0), ( "NO", 1.0),
                               ( "O2", 1.0), ( "C2", 1.0), ("NCO", 0.0), ("NO+", 0.0), ("O2+", 0.0),
			       ("CO+", 0.0), ("N2+", 0.0), ("CN+", 0.0), ( "Ar", 1.0), (  "C",  eA),
			       (  "N",  eA), (  "O",  eA), ( "C+", 0.0), ( "O+", 0.0), ( "e-", eE)  ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[2]  )

# 	1.3 O2 Dissociation reactions

eA = 1.000e22 / 2.000e21

r3 = make_reaction("O2 + M <=> O + O + M",
		fr=Park_NE_forward( 2.000e21, -1.5, 59750.0, 5, 0.5, 0.0 ),
                efficiencies=[ ("CO2", 1.0), ( "CO", 1.0), ( "N2", 1.0), ( "CN", 1.0), ( "NO", 1.0),
                               ( "O2", 1.0), ( "C2", 1.0), ("NCO", 0.0), ("NO+", 0.0), ("O2+", 0.0),
			       ("CO+", 0.0), ("N2+", 0.0), ("CN+", 0.0), ( "Ar", 1.0), (  "C",  eA),
			       (  "N",  eA), (  "O",  eA), ( "C+", 0.0), ( "O+", 0.0), ( "e-", 1.0)  ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[5]  )

# 	1.4 CN Dissociation reaction

C_CN = 2.530e14
    
r4 = make_reaction("CN + M <=> C + N + M",
		fr=Park_NE_forward( C_CN, 0.0, 71000.0, 3, 0.5, 0.0 ),
                reac_type="ThirdBodyReaction",
		participating_ivib=[3] )

# 	1.5 CO Dissociation reactions

eA  = 3.400e20 / 2.300e20
eAr = 2.300e19 / 2.300e20

r5 = make_reaction("CO + M <=> C + O + M",
		fr=Park_NE_forward( 2.300e20, -1.0, 129000.0, 1, 0.5, 0.0 ),
                efficiencies=[ ("CO2", 1.0), ( "CO", 1.0), ( "N2", 1.0), ( "CN", 1.0), ( "NO", 1.0),
                               ( "O2", 1.0), ( "C2", 1.0), ("NCO", 0.0), ("NO+", 0.0), ("O2+", 0.0),
			       ("CO+", 0.0), ("N2+", 0.0), ("CN+", 0.0), ( "Ar", eAr), (  "C",  eA),
			       (  "N",  eA), (  "O",  eA), ( "C+", 0.0), ( "O+", 0.0), ( "e-", 1.0)  ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[1]  )

# 	1.6 NO Dissociation reactions

C_NO  = 9.640e14
Ed_NO = 74700.0
eN2   = 1.450e15 / C_NO
eCO2  = 2.410e15 / C_NO
eX    = 1.0

r6 = make_reaction("NO + M <=> N + O + M",
		fr=Park_NE_forward( C_NO, 0.0, Ed_NO, 4, 0.5, 0.0 ),
                efficiencies=[ ("CO2",eCO2), ( "CO", 1.0), ( "N2", eN2), ( "CN", 1.0), ( "NO",  eX),
                               ( "O2", 1.0), ( "C2", 1.0), ("NCO", 0.0), ("NO+", 0.0), ("O2+", 0.0),
			       ("CO+", 0.0), ("N2+", 0.0), ("CN+", 0.0), ( "Ar", 1.0), (  "C",  eX),
			       (  "N",  eX), (  "O",  eX), ( "C+", 0.0), ( "O+", 0.0), ( "e-", 1.0)  ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[4]  )

# 	1.7 CO2 Dissociation reactions

C_CO2 = 6.900e21
eAr = 6.900e20 / C_CO2
eA  = 1.400e22 / C_CO2

r7 = make_reaction("CO2 + M <=> CO + O + M",
		fr=Park_NE_forward( C_CO2, -1.5, 63275.0, 0, 0.5, 0.0 ),
                efficiencies=[ ("CO2", 1.0), ( "CO", 1.0), ( "N2", 1.0), ( "CN", 1.0), ( "NO", 1.0),
                               ( "O2", 1.0), ( "C2", 1.0), ("NCO", 0.0), ("NO+", 0.0), ("O2+", 0.0),
			       ("CO+", 0.0), ("N2+", 0.0), ("CN+", 0.0), ( "Ar", eAr), (  "C",  eA),
			       (  "N",  eA), (  "O",  eA), ( "C+", 0.0), ( "O+", 0.0), ( "e-", 1.0)  ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[0]  )

#       1.8 NCO Dissociation reaction

r8 = make_reaction("NCO + M <=> CO + N + M",
    		fr=Park_NE_forward( 6.300e16, -0.5, 24000.0, 7, 0.5, 0.0 ),
		reac_type="ThirdBodyReaction",
		participating_ivib=[7] )

# 2.0 Exchange Reactions

r9 = make_reaction("N2 + O <=> NO + N",
		fr=GA_forward( 6.400e17, -1.0, 38370.0 ),
                reac_type="Exchange" )

r10 = make_reaction("NO + O <=> N + O2",
		fr=GA_forward( 8.400e12, 0.0, 19450.0 ),
                reac_type="Exchange" )

r11 = make_reaction("CO + C <=> C2 + O",
		fr=GA_forward( 2.000e17, -1.0, 58000.0 ),
                reac_type="Exchange" )

r12 = make_reaction("CO + O <=> C + O2",
		fr=GA_forward( 3.900e13, -0.18, 69200.0 ),
                reac_type="Exchange" )

r13 = make_reaction("CO + N <=> CN + O",
		fr=GA_forward( 1.000e14, 0.0, 38600.0 ),
                reac_type="Exchange" )

r14 = make_reaction("N2 + C <=> CN + N",
		fr=GA_forward( 5.240e13, 0.0, 22600.0 ),
                reac_type="Exchange" )

r15 = make_reaction("CN + O <=> NO + C",
		fr=GA_forward( 1.600e13, 0.10, 14600.0 ),
                reac_type="Exchange" )

r16 = make_reaction("CN + C <=> C2 + N",
		fr=GA_forward( 5.000e13, 0.0, 13000.0 ),
                reac_type="Exchange" )

# r17 is an addition reaction from Losev
r17 = make_reaction("CO + CO <=> C + CO2",
		fr=GA_forward( 2.300e9, 0.50, 65710.0 ),
                reac_type="Exchange" )

r18 = make_reaction("CO2 + O <=> O2 + CO",
		fr=GA_forward( 2.100e13, 0.0, 27800.0 ),
                reac_type="Exchange" )

# r19 is an additional reaction from Gokcen
r19 = make_reaction("C2 + N2 <=> CN + CN",
		fr=GA_forward( 1.500e13, 0.0, 21000.0 ),
                reac_type="Exchange" )

r20 = make_reaction("CO + NO <=> NCO + O",
		fr=GA_forward( 3.800e17, -0.873, 51600.0 ),
                reac_type="Exchange" )

r21 = make_reaction("CN + O2 <=> NCO + O",
		fr=GA_forward( 6.600e12,   0.00, -200.0  ),
                reac_type="Exchange" )

r22 = make_reaction("CN + CO2 <=> NCO + CO",
		fr=GA_forward( 4.000e14,   0.00, 19200.0 ),
                reac_type="Exchange" )

r23 = make_reaction("CN + NO <=> NCO + N",
		fr=GA_forward( 1.000e14,   0.00, 21200.0 ),
                reac_type="Exchange" )

r24 = make_reaction("CN + CO <=> NCO + C",
		fr=GA_forward( 1.500e16,  -4.87, 65800.0 ),
                reac_type="Exchange" )

# r25 is an addition reaction from Losev
r25 = make_reaction("CO + N <=> NO + C",
		fr=GA_forward( 2.900e11, 0.50, 53630.0 ),
                reac_type="Exchange" )

# r26 is an addition reaction from Losev
r26 = make_reaction("NO + CO <=> CO2 + N",
		fr=GA_forward( 4.600e8, 0.50, 12070.0 ),
                reac_type="Exchange" )

# 3.0 Ionic reactions

# 	3.1 Associatiove ionization reactions
r27 = make_reaction("N + O <=> NO+ + e-",
		fr = GA_forward( 8.8e8, 1.0, 31900.0 ) )

# 	3.2 Charge exchange reactions
r32 = make_reaction("NO+ + C <=> NO + C+",
		fr = GA_forward( 1.0e13, 0.0, 23200.0 ) )

r44 = make_reaction("C+ + e- <=> C",
    		fr = Park_NE_forward( 2.02e11, -0.46, 0.0, 0, 0.0, 1.0 ) )

# The old set of park reactions I used to use
reactions_list = [ r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9, r10,
                  r11, r12, r13, r14, r15, r16, r17, r18, r19, r20,
		  r21, r22, r23, r24, r25, r26, r27, r32, r44 ]

ndata.declare_reactions( reactions_list )
