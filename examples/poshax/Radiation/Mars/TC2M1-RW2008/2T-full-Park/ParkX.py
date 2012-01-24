# ParkX-NE-B.py
#
# Complete reaction scheme

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

# some useful flags
INCLUDE_LEE_MODS = False
PARK_NO_RATE = True
PARK_C2_RATE = True
INCLUDE_NCO = True
INCLUDE_PARK_EXTRAS = True
INCLUDE_GOKCEN_REACTIONS = True
INCLUDE_FUJITA_MODS = False
INCLUDE_ROND_MODS = False
INCLUDE_GORELOV_MODS = False
NO_RRRR = True

SEPARATE_EID_REACTIONS = True
if SEPARATE_EID_REACTIONS:
    fe = 0.0
else:
    fe = 1.0

ndata.declare_species( species_list )

# 1.0 Dissociation Reactions

# 	1.1 C2 Dissociation reaction

if INCLUDE_LEE_MODS==True and PARK_C2_RATE==False:
    # Lee's faster reaction rate [5x rate is TOO fast]
    # C_C2 = 7.500e16
    C_C2 = 1.500e16
else:
    C_C2 = 3.700e14

r1 = make_reaction("C2 + M <=> C + C + M",
    		fr=Park_NE_forward( C_C2, 0.0, 69900.0, 6, 0.5, 0.0 ),
		efficiencies= [ ("e-", fe) ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[6] )

r1e = make_reaction("C2 + e- <=> C + C + e-",
    		fr=Park_NE_forward( C_C2, 0.0, 69900.0, 6, 0.5, 0.5 ),
                reac_type="Dissociation",
		participating_ivib=[6] )

# 	1.2 N2 Dissociation reactions

eA = 3.000e22 / 7.000e21
if INCLUDE_LEE_MODS==True:
    # Gokcen's slower efficiency for electrons
    eE = 3.000e24 / 7.000e21
else:
    eE = 1.200e25 / 7.000e21

r2 = make_reaction("N2 + M <=> N + N + M",
		fr=Park_NE_forward( 7.000e21, -1.6, 113200.0, 2, 0.5, 0.0 ),
                efficiencies=[ ("CO2", 1.0), ( "CO", 1.0), ( "N2", 1.0), ( "CN", 1.0), ( "NO", 1.0),
                               ( "O2", 1.0), ( "C2", 1.0), ("NCO", 0.0), ("NO+", 0.0), ("O2+", 0.0),
			       ("CO+", 0.0), ("N2+", 0.0), ("CN+", 0.0), ( "Ar", 1.0), (  "C",  eA),
			       (  "N",  eA), (  "O",  eA), ( "C+", 0.0), ( "O+", 0.0), ( "e-", fe*eE)  ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[2]  )

r2e = make_reaction("N2 + e- <=> N + N + e-",
		fr=Park_NE_forward( eE*7.000e21, -1.6, 113200.0, 2, 0.5, 0.5 ),
                reac_type="Dissociation",
		participating_ivib=[2]  )

# 	1.3 O2 Dissociation reactions

eA = 1.000e22 / 2.000e21

r3 = make_reaction("O2 + M <=> O + O + M",
		fr=Park_NE_forward( 2.000e21, -1.5, 59750.0, 5, 0.5, 0.0 ),
                efficiencies=[ ("CO2", 1.0), ( "CO", 1.0), ( "N2", 1.0), ( "CN", 1.0), ( "NO", 1.0),
                               ( "O2", 1.0), ( "C2", 1.0), ("NCO", 0.0), ("NO+", 0.0), ("O2+", 0.0),
			       ("CO+", 0.0), ("N2+", 0.0), ("CN+", 0.0), ( "Ar", 1.0), (  "C",  eA),
			       (  "N",  eA), (  "O",  eA), ( "C+", 0.0), ( "O+", 0.0), ( "e-", 0.0)  ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[5]  )

# 	1.4 CN Dissociation reaction

if INCLUDE_LEE_MODS==True:
    # Gokcen's slightly faster reaction rate
    C_CN = 2.530e14
else:
    C_CN = 2.500e14
    
r4 = make_reaction("CN + M <=> C + N + M",
		fr=Park_NE_forward( C_CN, 0.0, 71000.0, 3, 0.5, 0.0 ),
		efficiencies= [ ("e-", fe) ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[3] )

r4e = make_reaction("CN + e- <=> C + N + e-",
		fr=Park_NE_forward( C_CN, 0.0, 71000.0, 3, 0.5, 0.5 ),
                reac_type="Dissociation",
		participating_ivib=[3] )

# 	1.5 CO Dissociation reactions

C_CO = 2.300e20
if INCLUDE_FUJITA_MODS:
    eA  = 2.000e20 / C_CO
    eAr = eA
    nCO = -1.2
else:
    eA  = 3.400e20 / C_CO
    eAr = 2.300e19 / C_CO
    nCO = -1.0

r5 = make_reaction("CO + M <=> C + O + M",
		fr=Park_NE_forward( C_CO, nCO, 129000.0, 1, 0.5, 0.0 ),
                efficiencies=[ ("CO2", 1.0), ( "CO", 1.0), ( "N2", 1.0), ( "CN", 1.0), ( "NO", 1.0),
                               ( "O2", 1.0), ( "C2", 1.0), ("NCO", 0.0), ("NO+", 0.0), ("O2+", 0.0),
			       ("CO+", 0.0), ("N2+", 0.0), ("CN+", 0.0), ( "Ar", eAr), (  "C",  eA),
			       (  "N",  eA), (  "O",  eA), ( "C+", 0.0), ( "O+", 0.0), ( "e-", 0.0)  ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[1]  )

# 	1.6 NO Dissociation reactions

if INCLUDE_LEE_MODS==True and PARK_NO_RATE==False:
    # Tsang's slightly slower values
    C_NO  = 9.640e14
    Ed_NO = 74700.0
    eN2   = 1.450e15 / C_NO
    eCO2  = 2.410e15 / C_NO
    eX    = 1.0
else:
    C_NO  = 5.000e15
    Ed_NO = 75500
    eN2   = 1.0
    eX    = 1.100e17 / C_NO
    eCO2  = eX

r6 = make_reaction("NO + M <=> N + O + M",
		fr=Park_NE_forward( C_NO, 0.0, Ed_NO, 4, 0.5, 0.0 ),
                efficiencies=[ ("CO2",eCO2), ( "CO", 1.0), ( "N2", eN2), ( "CN", 1.0), ( "NO",  eX),
                               ( "O2", 1.0), ( "C2", 1.0), ("NCO", 0.0), ("NO+", 0.0), ("O2+", 0.0),
			       ("CO+", 0.0), ("N2+", 0.0), ("CN+", 0.0), ( "Ar", 1.0), (  "C",  eX),
			       (  "N",  eX), (  "O",  eX), ( "C+", 0.0), ( "O+", 0.0), ( "e-", 0.0)  ],
                reac_type="ThirdBodyReaction",
		participating_ivib=[4]  )

# 	1.7 CO2 Dissociation reactions

if INCLUDE_ROND_MODS==True:
    # Blanket rate for all M's
    C_CO2 = 1.900e16
    eAr = 1.0
    eA  = 1.0
else:
    # standard Park (1994) rate(s)
    C_CO2 = 6.900e21
    eAr = 6.900e20 / C_CO2
    eA  = 1.400e22 / C_CO2

r7 = make_reaction("CO2 + M <=> CO + O + M",
		fr=Park_NE_forward( C_CO2, -1.5, 63275.0, 0, 0.5, 0.0 ),
                efficiencies=[ ("CO2", 1.0), ( "CO", 1.0), ( "N2", 1.0), ( "CN", 1.0), ( "NO", 1.0),
                               ( "O2", 1.0), ( "C2", 1.0), ("NCO", 0.0), ("NO+", 0.0), ("O2+", 0.0),
			       ("CO+", 0.0), ("N2+", 0.0), ("CN+", 0.0), ( "Ar", eAr), (  "C",  eA),
			       (  "N",  eA), (  "O",  eA), ( "C+", 0.0), ( "O+", 0.0), ( "e-", 0.0)  ],
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

if INCLUDE_FUJITA_MODS==True:
    r12 = make_reaction("CO + O <=> C + O2",
			fr=GA_forward( 1.000e15, -0.2, 70800.0 ),
			reac_type="Exchange" )	
else:
    r12 = make_reaction("CO + O <=> C + O2",
			fr=GA_forward( 3.900e13, -0.18, 69200.0 ),
			reac_type="Exchange" )

if INCLUDE_ROND_MODS==True:
    r13 = make_reaction("CO + N <=> CN + O",
		fr=GA_forward( 1.000e12, -0.5, 38600.0 ),
                reac_type="Exchange" )
else:
    r13 = make_reaction("CO + N <=> CN + O",
		fr=GA_forward( 1.000e14, 0.0, 38600.0 ),
                reac_type="Exchange" )

if INCLUDE_LEE_MODS==True:
    r14 = make_reaction("N2 + C <=> CN + N",
		fr=GA_forward( 5.240e13, 0.0, 22600.0 ),
                reac_type="Exchange" )
else:
    r14 = make_reaction("N2 + C <=> CN + N",
		fr=GA_forward( 1.100e14, -0.11, 23200.0 ),
                reac_type="Exchange" )

if INCLUDE_ROND_MODS==True:
    r15 = make_reaction("CN + O <=> NO + C",
		fr=GA_forward( 1.600e11, 0.10, 14600.0 ),
                reac_type="Exchange" )
else:
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

# Unsure about the implementation of this reaction...
r28 = make_reaction("O + O <=> O2+ + e-",
		fr = GA_forward( 7.1e2, 2.7, 80600.0 ) )

r29 = make_reaction("C + O <=> CO+ + e-",
		fr = GA_forward( 8.8e8, 1.0, 33100.0 ) )

# N2+ forming reaction from Gokcen
r30 = make_reaction("N + N <=> N2+ + e-",
                fr=GA_forward( 4.40e7, 1.50, 67500.0 ) )

# CN+ forming reaction from Gokcen [unstable]
r31 = make_reaction("C + N <=> CN+ + e-",
                fr=GA_forward( 1.00e15, 1.50, 164400.0 ) )

# 	3.2 Charge exchange reactions
r32 = make_reaction("NO+ + C <=> NO + C+",
		fr = GA_forward( 1.0e13, 0.0, 23200.0 ) )

r33 = make_reaction("O2+ + O <=> O+ + O2",
		fr = GA_forward( 4.0e12, -0.09, 18000.0 ) )

r34 = make_reaction("NO+ + N <=> O+ + N2",
		fr = GA_forward( 3.4e13, -1.08, 12800.0 ) )

r35 = make_reaction("NO+ + O <=> O2+ + N",
		fr = GA_forward( 7.2e12, 0.29, 48600.0 ) )

if INCLUDE_GORELOV_MODS==True:
    r36 = make_reaction("CO + C+ <=> CO+ + C",
		fr = GA_forward( 1.0e12, 0.0, 31400.0 ) )
else:
    # Park 1994 rate
    r36 = make_reaction("CO + C+ <=> CO+ + C",
		fr = GA_forward( 1.0e13, 0.0, 31400.0 ) )

if INCLUDE_GORELOV_MODS==True:
    r37 = make_reaction("O2 + C+ <=> O2+ + C",
		fr = GA_forward( 1.0e12, 0.0, 9400.0 ) )
else:
    # Park 1994 rate
    r37 = make_reaction("O2 + C+ <=> O2+ + C",
		fr = GA_forward( 1.0e13, 0.0, 9400.0 ) )

# N2+ forming reaction from Gokcen
r38 = make_reaction("C+ + N2 <=> N2+ + C",
                fr=GA_forward( 1.11e14, -0.11, 50000.0 ) )

# CN+ consuming reaction from Gokcen
r39 = make_reaction("CN+ + N <=> CN + N+",
                fr=GA_forward( 9.80e12, 0.0, 40700.0 ) )

# 	3.3 Electron-impact ionization
#           NOTE: Using T_e as rate controlling temperature

r40 = make_reaction("C + e- <=> C+ + e- + e-",
		fr = Park_NE_forward( 3.9e33, -3.78, 130700.0, 0, 0.0, 1.0 ) )

# N+ forming reaction from Gokcen
r41 = make_reaction("N + e- <=> N+ + e- + e-",
    		fr = Park_NE_forward( 2.500e34, -3.82, 168600.0, 0, 0.0, 1.0 ) )

r42 = make_reaction("O + e- <=> O+ + e- + e-",
		fr = Park_NE_forward( 3.9e33, -3.78, 158500.0, 0, 0.0, 1.0 ) )

# 	3.4 Radiative recombination reactions
#           NOTE: Using T_e as rate controlling temperature

r43 = make_reaction("O+ + e- <=> O",
    		fr = Park_NE_forward( 1.07e11, -0.52, 0.0, 0, 0.0, 1.0 ) )

r44 = make_reaction("C+ + e- <=> C",
    		fr = Park_NE_forward( 2.02e11, -0.46, 0.0, 0, 0.0, 1.0 ) )
		
if NO_RRRR:
    r43.set_br_model( Park_NE_backward( 0.0, 0.0, 0.0, 0, 0.0, 1.0 ) )
    r44.set_br_model( Park_NE_backward( 0.0, 0.0, 0.0, 0, 0.0, 1.0 ) )

# The old set of park reactions I used to use
reactions_list = [ r1,  r2,  r3,  r4,  r5,  r6,  r7,  r9,
                  r10, r11, r12, r13, r14, r15, r16, r18, r27,
		  r29, r32, r33, r34, r35, r36, r37, r40 ]

# Park NCO reactions that I previously left out
NCO_reactions = [  r8, r20, r21, r22, r23, r24 ]

# additional reactions included in Lee's scheme
new_lee_reactions = [ r17, r19, r25, r26 ]

# O2+ association ionisation, rad recombination (previously omitted)
park_extras = [ r28, r42, r43, r44 ]

# Gokcen reactions: r31 is unstable
gokcen_reactions = [ r30, r38, r39, r41 ]

# Separately defined electron impact dissociation reactions
eid_reactions = [ r1e, r2e, r4e ]

if INCLUDE_NCO:
    reactions_list += NCO_reactions

if INCLUDE_LEE_MODS:
    reactions_list += new_lee_reactions

if INCLUDE_PARK_EXTRAS:
    reactions_list += park_extras
    
if INCLUDE_GOKCEN_REACTIONS:
    reactions_list += gokcen_reactions
    
if SEPARATE_EID_REACTIONS:
    reactions_list += eid_reactions

ndata.declare_reactions( reactions_list )
