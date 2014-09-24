import sys

class ParadeRadiator:
    """Base parade radiator class"""
    def __init__( self, name, atoms, mw, dat_file, systems=[] ):
        self.pname = name
	self.name = name.replace("+","_plus").replace("E","e_minus")
        self.atoms = atoms
	if self.atoms==0:
	    self.type = "electron_radiator"
	elif self.atoms==1:
	    self.type = "atomic_radiator"
	elif self.atoms==2:
	    self.type = "diatomic_radiator"
	elif self.atoms==3:
	    self.type = "triatomic_radiator"
	else:
	    print "atoms: %d not understood"
	    sys.exit()
	self.systems = systems
	self.dat_file = dat_file
	self.mol_weight = mw *1.0e-3
        self.band_model = 'LBL'
	self.isp = -1
	self.iT = 0
	self.iTe = 0
	self.iTv = 0
	self.iTr = 0

    def get_parade_string( self, mol_index, data_path ):
        if data_path[-1]!="/":
            data_path += "/"
        data_file_path = data_path + self.dat_file
        #if self.atoms>1:
        #    ttr = mol_index*2 + 1
        #    tv  = mol_index*2 + 2
        #else:
        #    ttr = 1
        #    tv  = 2  
	## this has been updated based on the fireii example in PARADE v3.1.
        if self.atoms==0:
	    tt = self.iTe+1
	    tr = 0
            tv = 0
	elif self.atoms==1:
	    tt = self.iT+1
            tr = 0
            tv  = 0
        else:
            tt = self.iT+1
            tr = self.iTr+1
            tv  = self.iTv+1
	atoms = self.atoms
        if self.band_model == "SNB":
            atoms *= -1              
        string = " 'Y',        '%s'           %d       %d   %d   %d   2,   '%s'\n" % ( self.pname, atoms, tt, tr, tv, data_file_path )
        for sys in self.systems:
            string += " 'b'    '%s' \n" % sys		# sys is the string for the band system. random 'N' was removed.
        return string

    def get_LUA_string(self):
	string  = ""
	string += "%s = {}\n" % ( self.name )
	string += "%s.isp = %d\n" % ( self.name, self.isp )
	string += "%s.type = '%s'\n" % ( self.name, self.type )
	string += "%s.mol_weight = %e\n" % ( self.name, self.mol_weight )
	if self.type=="diatomic_radiator" or self.type=="triatomic_radiator":
	    string += "%s.iTv = %d\n" % ( self.name, self.iTv )
	    string += "%s.iTr = %d\n" % ( self.name, self.iTr )
	return string
        
    def default_data( self ):
        # nothing to do for the moment
        return

# when checking the species data files, check also the specified radiating bands.
available_radiators_boltzmann = { "Ar"       :    ParadeRadiator( "Ar" , 1, 39.948   ,    "NIST/arinist.dat" ),
                        "Ar_plus"  :    ParadeRadiator( "Ar+", 1, 39.948   ,   "NIST/ariinist.dat" ),
                        "C"        :    ParadeRadiator( "C"  , 1, 12.011   ,     "NIST/cinistmf.dat" ),
                        "F"        :    ParadeRadiator( "F"  , 1, 18.9984  ,     "NIST/finist.dat" ),
                        "H"        :    ParadeRadiator( "H"  , 1,  1.008   ,     "NIST/hinist.dat" ),
                        "H_plus"   :    ParadeRadiator( "H+" , 1,  1.008   ,    "hii.dat" ),
                        "He"       :    ParadeRadiator( "He" , 1,  4.0026  ,    "NIST/heinist.dat" ),
                        "N"        :    ParadeRadiator( "N"  , 1, 14.0067  ,     "NIST/niNISTmf.dat" ),
                        "N_plus"   :    ParadeRadiator( "N+" , 1, 14.0067  ,    "NIST/niiNISTmf.dat" ),
                        "O"        :    ParadeRadiator( "O"  , 1, 15.998   ,     "NIST/oiNISTmf.dat" ),
                        "O_plus"   :    ParadeRadiator( "O+" , 1, 15.998   ,    "oiiseb.dat" ),		# changed from oiiNISTmf.dat
                        "C2"       :    ParadeRadiator( "C2" , 2, 24.0214  ,  "C2ijb.dat", [ "swan" ] ),
                        "CH"       :    ParadeRadiator( "CH" , 2, 13.0186  ,  "CHibp.dat", [ "3900A", "4300A" ] ),
                        "CN"       :    ParadeRadiator( "CN" , 2, 26.0174  , "CNibp.dat", [ "CNviolet", "CNred" ] ),	# may not have CN red... check!
                        "CO"       :    ParadeRadiator( "CO" , 2, 28.0101  ,  "COilh.dat", [ "CO4p", "CO3p", "COA", "Asundi" ] ),
                        "CO_plus"  :    ParadeRadiator( "CO+", 2, 28.0101  , "COiimf.dat", [ "COp1n", "COpComet" ] ),
                        "N2"       :    ParadeRadiator( "N2" , 2, 28.01348 ,  "N2ihl.dat", [ "N21p", "N22p", "N2bh2", "N2bh", "N2cy", "N2wj", "N2w", "N2eX" ] ),
                        "N2_plus"  :    ParadeRadiator( "N2+", 2, 28.012931, "N2iibp.dat", [ "N2p1n" ] ),
                        "NH"       :    ParadeRadiator( "NH" , 2, 15.01468 , "NHibp.dat" , [ "3360A" ] ),
                        "NO"       :    ParadeRadiator( "NO" , 2, 30.00614 , "NOibp.dat" , [ "NObeta", "NOgamma", "NOdelta", "NOepsilon" ] ),
                        "O2"       :    ParadeRadiator( "O2" , 2, 31.9988  , "O2ibp.dat" , [ "O2sr" ] ),
			"H2"	   :    ParadeRadiator( "H2" , 2, 2.016    , "H2inj.dat" ),
                        "CO2"      :    ParadeRadiator( "CO2", 3, 44.0095  , "CO2ihl.dat" , [ "CO2rovib" ] ),
			"C3"	   :    ParadeRadiator( "C3" , 3, 36.033   , "C3_FGE.cs" ),
			"C2H"      :    ParadeRadiator( "C2H", 3, 25.030   , "C2H_KAIST.cs" ),
                        "e_minus"  :    ParadeRadiator( "E" , 0, 0.000548579903, "none" ) }

# remember to include important bands, like the N2 Rydberg bands...
available_radiators_qss = { "Ar"       :    ParadeRadiator( "Ar" , 1, 39.948   , "NIST/arinist.dat"  ),
						"Ar_plus"  :    ParadeRadiator( "Ar+", 1, 39.948   ,  "NIST/ariinist.dat" ),
						"C"        :    ParadeRadiator( "C"  , 1, 12.011   ,  "NIST/cinistmf.dat" ),
						"F"        :    ParadeRadiator( "F"  , 1, 18.9984  ,  "NIST/finist.dat" ),
						"H"        :    ParadeRadiator( "H"  , 1,  1.008   ,  "NIST/hinist.dat" ),
						"H_plus"   :    ParadeRadiator( "H+" , 1,  1.008   ,  "hii.dat" ),
						"He"       :    ParadeRadiator( "He" , 1,  4.0026  ,  "NIST/heinist.dat" ),
						"N"        :    ParadeRadiator( "N"  , 1, 14.0067  ,  "QSS/NiParknist.qss" ),
						"N_plus"   :    ParadeRadiator( "N+" , 1, 14.0067  ,  "NIST/niiNISTmf.dat" ),
						"O"        :    ParadeRadiator( "O"  , 1, 15.998   ,  "QSS/OiPark.qss" ),
						"O_plus"   :    ParadeRadiator( "O+" , 1, 15.998   ,  "oiiseb.dat" ),
						"C2"       :    ParadeRadiator( "C2" , 2, 24.0214  ,  "C2ijb.dat", [ "swan" ] ),
						"CH"       :    ParadeRadiator( "CH" , 2, 13.0186  ,  "CHibp.dat", [ "3900A", "4300A" ] ),
						"CN"       :    ParadeRadiator( "CN" , 2, 26.0174  ,  "QSS/CNijb.qss", [ "CNviolet" ] ),
						"CO"       :    ParadeRadiator( "CO" , 2, 28.0101  ,  "COilh.dat", [ "CO4p", "CO3p", "COA", "Asundi" ] ),
						"CO_plus"  :    ParadeRadiator( "CO+", 2, 28.0101  ,  "COiimf.dat", [ "COp1n", "COpComet" ] ),
						"N2"       :    ParadeRadiator( "N2" , 2, 28.01348 ,  "QSS/N2inj.qss" , [ "N21p", "N22p", "N2bh2", "N2bh", "N2cy", "N2wj", "N2w", "N2eX" ] ),
						"N2_plus"  :    ParadeRadiator( "N2+", 2, 28.012931,  "QSS/N2iijb.qss", [ "N2p1n" ] ),
						"NH"       :    ParadeRadiator( "NH" , 2, 15.01468 ,  "NHibp.dat" , [ "3360A" ] ),
						"NO"       :    ParadeRadiator( "NO" , 2, 30.00614 ,  "QSS/NOijb.qss" , [ "NObeta" , "NOgamma" , "NOdelta" , "NOepsilon" ] ),
						"O2"       :    ParadeRadiator( "O2" , 2, 31.9988  ,  "QSS/O2ijb.qss" , [ "O2sr" ] ),
						"H2"	   :    ParadeRadiator( "H2" , 2, 2.016    ,  "H2inj.dat" ),
			                        "CO2"      :    ParadeRadiator( "CO2", 3, 44.0095  ,  "CO2ihl.dat" , [ "CO2rovib" ] ),
						"C3"	   :    ParadeRadiator( "C3" , 3, 36.033   ,  "C3_FGE.cs" ),
						"C2H"      :    ParadeRadiator( "C2H", 3, 25.030   ,  "C2H_KAIST.cs" ),
						"e_minus"  :    ParadeRadiator( "E" , 0, 0.000548579903, "none" ) }
