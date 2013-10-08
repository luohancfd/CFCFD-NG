import os
import sys

def add_commas( string ):
    tks = string.split()
    for tk in tks[:-1]:
        string = string.replace(tk+" ",tk+", ")
    return string
    
tab = "   "

class AtomicLevel(object):
    """Class describing an individual or multiplet atomic level"""
    def __init__(self,n=-1,E=-1,g=-1,l=-1,L=-1,S=-1,parity=-1,conf="",term=""):
        self.n = n
        self.E = E
        self.g = g
        self.l = l
        self.L = L
        self.S = S
        self.parity = parity
        self.conf = conf
        self.term = term
        
    def get_string(self):
        lev_string = "n = " + str(self.n)
        lev_string += ", E = " + str(self.E)
        lev_string += ", g = " + str(self.g)
        lev_string += ", l = " + str(self.l)
        lev_string += ", L = " + str(self.L)
        lev_string += ", S = " + str(self.S)
        lev_string += ", parity = " + str(self.parity)
        lev_string += ", conf = " + str(self.conf)
        lev_string += ", term = " + str(self.term)
        return lev_string

    def get_LUA_string(self):
        lua_lev_string = "%d, %e, %d, %d, %d, %d, %d" % ( self.n, self.E, \
                                   self.g, self.l, self.L, self.S, self.parity )
        return lua_lev_string

class AtomicLine:
    def __init__(self, g_u, E_u, g_l, E_l, A_ul, lambda_ul, conf_u="", term_u="", conf_l="", term_l="", acc="", type_str="", ilev_u=-1, ilev_l=-1, n=-1, gamma_S0=-1 ):
        self.g_u = g_u
        self.E_u = E_u
        self.g_l = g_l
        self.E_l = E_l
        self.A_ul = A_ul
        if lambda_ul > 0.0:
            self.lambda_ul = lambda_ul
        else:
            self.lambda_ul = 1.0e7 / ( self.E_u - self.E_l )
        self.conf_u = conf_u
        self.term_u = term_u
        self.conf_l = conf_l
        self.term_l = term_l
        self.acc = acc
        self.type_str = type_str
        # FIXME: make type from type_str
        self.type = -1
        self.ilev_u = ilev_u
        self.ilev_l = ilev_l
        self.n  = n
        self.gamma_S0 = gamma_S0
        
    def get_string(self):
        line_string = "g_u = " + str(self.g_u)
        line_string += ", E_u = " + str(self.E_u)
        line_string += ", g_l = " + str(self.g_l)
        line_string += ", E_l = " + str(self.E_l)
        line_string += ", A_ul = " + str(self.A_ul)
        line_string += ", lambda_ul = " + str(self.lambda_ul)
        line_string += ", conf_u = " + str(self.conf_u)
        line_string += ", term_u = " + str(self.term_u)
        line_string += ", conf_l = " + str(self.conf_l)
        line_string += ", term_l = " + str(self.term_l)
        line_string += ", acc = " + str(self.acc)
        line_string += ", type = " + str(self.type)
        return line_string

#----------------------------------------------------------------------------------
#                No.        Ei(cm-1)   Ek(cm-1)   gi   gk Aki(1/s)  ie_i ie_k type
#----------------------------------------------------------------------------------

    def get_LUA_string(self):
        lua_line_string = "%e, %e, %d, %d, %e, %d, %d, %d, %e, %e" % ( self.E_l, self.E_u, \
                                   self.g_l, self.g_u, self.A_ul, self.ilev_l, \
                                   self.ilev_u, self.type, self.n, self.gamma_S0 )
        return lua_line_string

class PhotoIonXSectionModel(object):
    """Base class for describing a photo-ionization cross-section model"""
    def __init__(self, elevel_set="", model=""):
            self.model = model
            self.elevel_set = elevel_set
            self.comments = "# Description of the photoionization cross-section model"
            
    def get_LUA_string(self, species):
        ostring  = ""
        ostring += "%s.photoionXsection_model = {\n" % ( species )
        comments = self.comments.replace('#','   --')
        ostring += "%s\n" % ( comments )
        ostring += tab+"model = '%s',\n" % ( self.model )
        ostring += "}\n"
        return ostring
        
class JohnstonModel(PhotoIonXSectionModel):
    """Derived class for describing Johnston's photo-ionization cross-section model"""
    def __init__(self, elevel_set=""):
        PhotoIonXSectionModel.__init__(self, elevel_set=elevel_set, model="JohnstonModel")
        self.steps = []
        self.thresholds = []
            
    def get_LUA_string(self, species):
        ostring  = ""
        ostring += "%s.photoionXsection_model = {\n" % ( species )
        comments = self.comments.replace('#','   --')
        ostring += "%s\n" % ( comments )
        ostring += tab+"model = '%s',\n" % ( self.model )
        ostring += tab+"nsteps = %d,\n" % ( len(self.steps) )
        tks = self.steps[-1].split()
        nstep_levs = int(tks[0])+1
        ostring += tab+"nstep_levs = %d,\n" % ( nstep_levs )
        ostring += tab+"-- ====================================================================\n"
        ostring += tab+"--   No.     ilev    E_a(eV)    E_b(eV)      sigma_bf x 10^18 (cm^2)   \n"
        ostring += tab+"-- ====================================================================\n"
        for istep,step in enumerate(self.steps):
            if istep<10: ostring += tab+"step_%d   = { %s },\n" % (istep,add_commas(step))
            else: ostring += tab+"step_%d  = { %s },\n" % (istep,add_commas(step))
        ostring += tab+"-- ====================================================================\n"
        ostring += tab+"nthresholds = %d,\n" % ( len(self.thresholds) )
        ostring += tab+"-- =======================================================================\n"
        ostring += tab+"--   No.           ilev    E (eV)     sigma_bf x 10^18 (cm^2)     theta   \n"
        ostring += tab+"-- =======================================================================\n"
        for ithreshold,threshold in enumerate(self.thresholds):
            if ithreshold<10: ostring += tab+"threshold_%d   = { %s },\n" % ( ithreshold, add_commas(threshold) )
            else: ostring += tab+"threshold_%d  = { %s },\n" % ( ithreshold, add_commas(threshold) )
        ostring += tab+"-- =======================================================================\n"
        ostring += "}\n"
        return ostring

class TOPBasePICSLevel(object):
    """Derived class for describing TOPBase photo-ionization cross-section model"""
    def __init__(self, E, term, ilevTB, E_Ryd_list=[], sigma_list=[]):
        self.E = E
        self.term = term
        self.ilevTB = ilevTB
        self.E_Ryd_list = E_Ryd_list
        self.sigma_list = sigma_list

    def get_LUA_string(self):
        ostring = ""
        ostring += tab + "ilev_%d = {\n" % self.ilev
        ostring += tab + tab + "E = %e,\n" % self.E
        ostring += tab + tab + "npoints = %d,\n" % len(self.E_Ryd_list)
        ostring += tab +tab+"-- =======================================================================\n"
        ostring += tab +tab+"--   No.          E (Ry)     sigma_bf (cm^2)                   \n"
        ostring += tab +tab+"-- =======================================================================\n"
        for i in range(len(self.E_Ryd_list)):
            ostring += tab + tab + "point_%d = { %e, %e },\n" % ( i, self.E_Ryd_list[i], self.sigma_list[i] )
        ostring += tab + "},\n"
        return ostring
        
class TOPBasePICSModel(PhotoIonXSectionModel):
    """Derived class for describing TOPBase photo-ionization cross-section model"""
    def __init__(self, level_data, hydrogenic_fill=False, elevel_set="TOPBase"):
        PhotoIonXSectionModel.__init__(self, elevel_set=elevel_set, model="TOPBaseModel")
        self.level_data = level_data
        self.hydrogenic_fill = hydrogenic_fill

    def get_LUA_string(self, species):
        ostring  = ""
        ostring += "%s.photoionXsection_model = {\n" % ( species )
        comments = self.comments.replace('#','   --')
        ostring += "%s\n" % ( comments )
        ostring += tab+"model = '%s',\n" % ( self.model )
        ostring += tab + "hydrogenic_fill = %s,\n" % ( str(self.hydrogenic_fill).lower() )
        for ipl,level in enumerate(self.level_data):
            ostring += level.get_LUA_string()
        ostring += tab+"-- =======================================================================\n"
        ostring += "}\n"
        return ostring
        
class AtomicQSSModel(object):
    """Class for describing an atomic QSS model"""
    def __init__(self, name="", noneq_elevs="", eie_model="none", eii_model="none", rt_model="none", pr_model="none", special=""):
            self.name = name
            self.noneq_elevs = noneq_elevs
            self.eie_model = eie_model
            self.eii_model = eii_model
            self.rt_model = rt_model
            self.pr_model = pr_model
            self.special = special
            self.comments = "# Description of the atomic QSS model"
            self.inc_eq_elevs = 1
            self.T_lower = 4000.0
            
    def get_LUA_string(self, aname, special = ""):
        ostring  = ""
        ostring += "%s.QSS_model = {\n" % ( aname )
        comments = self.comments.replace('#','   --')
        ostring += "%s\n" % ( comments )
        ostring += tab+"noneq_elevs = { %s },\n" % ( self.noneq_elevs )
        ostring += tab+"inc_eq_elevs = %d,\n" % ( self.inc_eq_elevs )
        ostring += tab+"T_lower = %e,\n" % ( self.T_lower )
        ostring += tab+"electron_impact_excitation = '%s',\n" % ( self.eie_model )
        ostring += tab+"electron_impact_ionization = '%s',\n" % ( self.eii_model )
        ostring += tab+"radiative_transitions = '%s',\n" % ( self.rt_model )
        ostring += tab+"photorecombination = '%s',\n" % ( self.pr_model )
        ostring += "%s" % self.special
        ostring += "}\n"
        return ostring
        
class DiatomicQSSModel(object):
    """Class for describing a diatomic QSS model"""
    def __init__(self, name="", noneq_elevs="", noneq_elev_labels="", reactions=[]):
            self.name = name
            self.noneq_elevs = noneq_elevs
            self.noneq_elev_labels = noneq_elev_labels
            self.reactions = reactions
            self.comments = "# Description of the diatomic QSS model"
            self.inc_eq_elevs = 1
            self.T_lower = 4000.0
            
    def get_LUA_string(self, mname):
        ostring  = ""
        ostring += "%s.QSS_model = {\n" % ( mname )
        comments = self.comments.replace('#','   --')
        ostring += "%s\n" % ( comments )
        ostring += tab+"noneq_elevs = { %s },\n" % ( self.noneq_elevs )
        ostring += tab+"noneq_elev_labels = { %s },\n" % ( self.noneq_elev_labels )
        ostring += tab+"inc_eq_elevs = %d,\n" % ( self.inc_eq_elevs )
        ostring += tab+"T_lower = %e,\n" % ( self.T_lower )
        ostring += tab+"reactions = {\n"
        for reaction in self.reactions:
            ostring += tab*2+"{\n"
            for line in reaction:
                    ostring += tab*3 + line + "\n"
            ostring += tab*2+"},\n"
        ostring += tab + "},\n"
        ostring += "}\n"
        return ostring
        
class Radiator(object):
    """Base radiator class"""
    def __init__(self, name="", type=""):
        self.name = name
        self.type = type
        self.mol_weight = 0.0
        self.eta_I = 0.0
        self.h_f = 0.0
        self.Z = 0
        self.iT = 0
        self.iTe = 0
        self.E_pop_method = "boltzmann"
        self.isp = -1
        self.comments = "# Description of the radiator"
        
    def get_LUA_string(self):
        ostring  = ""
        ostring += "%s = {}\n" % ( self.name )
        comments = self.comments.replace('#','--')
        ostring += "%s\n" % ( comments )
        ostring += "%s.isp = %d\n" % ( self.name, self.isp )
        ostring += "%s.type = '%s'\n" % ( self.name, self.type )
        ostring += "%s.mol_weight = %e\n" % ( self.name, self.mol_weight )
        ostring += "%s.h_f = %e\n" % ( self.name, self.h_f )
        ostring += "%s.eta_I = %e\n" % ( self.name, self.eta_I )
        ostring += "%s.Z = %d\n" % ( self.name, self.Z )
        ostring += "%s.E_pop_method = '%s'\n" % ( self.name, self.E_pop_method )
        ostring += "%s.iT = %d\n" % ( self.name, self.iT )
        ostring += "%s.iTe = %d\n" % ( self.name, self.iTe )
        return ostring
        
class ElectronRadiator(Radiator):
    """Derived electron radiator class"""
    def __init__(self, name=""):
        Radiator.__init__(self, name, "electron_radiator")
        self.systems = []
        self.available_systems = {}
        
    def get_LUA_string(self):
        ostring  = ""
        ostring += "%s = {}\n" % ( self.name )
        comments = self.comments.replace('#','--')
        ostring += "%s\n" % ( comments )
        ostring += "%s.isp = %d\n" % ( self.name, self.isp )
        ostring += "%s.type = '%s'\n" % ( self.name, self.type )
        ostring += "%s.mol_weight = %e\n" % ( self.name, self.mol_weight )
        ostring += "%s.h_f = %e\n" % ( self.name, self.h_f )
        ostring += "%s.eta_I = %e\n" % ( self.name, self.eta_I )
        ostring += "%s.Z = %d\n" % ( self.name, self.Z )
        ostring += "%s.E_pop_method = 'none'\n" % self.name
        ostring += "%s.iT = %d\n" % ( self.name, self.iT )
        ostring += "%s.iTe = %d\n" % ( self.name, self.iTe )
        ostring += "%s.systems_list = { " % ( self.name )
        if len(self.systems)>0:
            for isys in self.systems:
                ostring += "'%s', " % ( isys )
        ostring += "}\n"
        return ostring
        
    def default_data( self ):
        for system in self.available_systems.keys():
             self.systems.append( self.available_systems[system] )
        
class AtomicRadiator(Radiator):
    """Derived atomic radiator class"""
    def __init__(self, name=""):
        Radiator.__init__(self, name, "atomic_radiator")
        self.line_set = AtomicLineSet()
        self.level_set = AtomicLevelSet()
        self.photoionXsection_model = PhotoIonXSectionModel()
        self.QSS_model = AtomicQSSModel()
        self.available_line_sets = {}
        self.available_level_sets = {}
        self.available_photoionXsection_models = {}
        self.available_QSS_models = {}
        self.default_line_set = ""
        self.default_level_set = ""
        self.default_photoionXsection_model = ""
        self.default_QSS_model = ""
        
    def get_LUA_string(self):
        ostring = Radiator.get_LUA_string(self)
        ostring += self.level_set.get_LUA_string(aname=self.name)
        ostring += self.line_set.get_LUA_string(aname=self.name)
        # FIXME: need to check that the level set is compatible with PICS model
        ostring += self.photoionXsection_model.get_LUA_string(species=self.name)
        # FIXME: need to check that the level set is compatible with the QSS model
        if self.E_pop_method=="QSS":
            ostring += self.QSS_model.get_LUA_string(aname=self.name)
        return ostring
        
    def default_data( self ):
        self.line_set = self.available_line_sets[self.default_line_set]
        self.level_set = self.available_level_sets[self.default_level_set]
        self.photoionXsection_model = self.available_photoionXsection_models[self.default_photoionXsection_model]
        self.QSS_model = self.available_QSS_models[self.default_QSS_model]
        
class AtomicLineSet(object):
    """Atomic line set class"""
    def __init__(self,lines=[],comments="# Description of the atomic line set",npoints=100,nwidths=1000,beta=1.01):
        self.lines = lines
        if comments[0]!="#": comments = "# " + comments
        self.comments = comments
        self.npoints = npoints
        self.nwidths = nwidths
        self.beta = beta
        
    def get_LUA_string(self, aname):
        ostring  = "%s.line_data = {\n" % ( aname )
        comments = self.comments.replace('#','   --')
        ostring += "%s\n" % ( comments )
        ostring += tab+"n_points = %d,\n" % ( self.npoints )
        ostring += tab+"n_widths = %d,\n" % ( self.nwidths )
        ostring += tab+"beta = %f,\n" % ( self.beta )
        ostring += tab+"n_lines = %d,\n" % ( len(self.lines) )
        ostring += tab+"-- ============================================================================\n"
        ostring += tab+"--    No.         Ei(cm-1)  Ek(cm-1)    gi    gk  Aki(1/s)    ie_i  ie_k  type \n"
        ostring += tab+"-- ============================================================================\n"
        for i,line in enumerate(self.lines):
            # FIXME: eventually all level data will be AtomicLevel instances and 
            #        these if and else's should be removed
            if type(line)==type(""):
                if i<10: ostring += tab+"iline_%d   = { %s },\n" % ( i, add_commas(line) )
                elif i<100: ostring += tab+"iline_%d  = { %s },\n" % ( i, add_commas(line) )
                else: ostring += tab+"iline_%d = { %s },\n" % ( i, add_commas(line) )
            else:
                if i<10: ostring += tab+"iline_%d   = { %s },\n" % ( i, line.get_LUA_string() )
                elif i<100: ostring += tab+"iline_%d  = { %s },\n" % ( i, line.get_LUA_string() )
                else: ostring += tab+"iline_%d = { %s },\n" % ( i, line.get_LUA_string() )
        ostring += tab+"-- ============================================================================\n"
        ostring += "}\n"
        return ostring
        
class AtomicLevelSet(object):
    """Atomic level set class"""
    def __init__(self,levels=[],comments="# Description of the atomic level set",isp_list=[]):
        self.levels = levels
        if comments[0]!="#": comments = "# " + comments
        self.comments = comments
        self.isp_list = isp_list
    def get_LUA_string(self, aname):
        ostring  = "%s.level_data = {\n" % ( aname )
        comments = self.comments.replace('#','   --')
        ostring += "%s\n" % ( comments )
        ostring += tab+"n_levels = %d,\n" % ( len(self.levels) )
        ostring += tab+"isp_list = { "
        for isp in self.isp_list:
            ostring += "%d, " % isp
        ostring += "},\n"
        ostring += tab+"-- ===========================================================\n"
        ostring += tab+"--   No.     n      E(cm-1)      g     l     L     S    parity \n"
        ostring += tab+"-- ===========================================================\n"
        for i,level in enumerate(self.levels):
            # FIXME: eventually all level data will be AtomicLevel instances and 
            #        these if and else's should be removed
            if type(level)==type(""):
                if i<10: ostring += tab+"ilev_%d  = { %s },\n" % ( i, add_commas(level) )
                else: ostring += tab+"ilev_%d = { %s },\n" % ( i, add_commas(level) )
            else:
                if i<10: ostring += tab+"ilev_%d  = { %s },\n" % ( i, level.get_LUA_string() )
                else: ostring += tab+"ilev_%d = { %s },\n" % ( i, level.get_LUA_string() )
        ostring += tab+"-- ===========================================================\n"
        ostring += "}\n"
        return ostring

class DiatomicRadiator(Radiator):
    """Derived diatomic radiator class"""
    def __init__(self, name=""):
        Radiator.__init__(self, name, "diatomic_radiator")
        self.systems = []
        self.level_set = DiatomicLevelSet()
        self.photoionXsection_model = PhotoIonXSectionModel()
        self.QSS_model = DiatomicQSSModel()
        self.available_systems = {}
        self.available_level_sets = {}
        self.available_photoionXsection_models = {}
        self.available_QSS_models = {}
        self.default_level_set = ""
        self.default_photoionXsection_model = ""
        self.default_QSS_model = ""
        self.iTv = 0
        self.iTr = 0
        self.I_spin = 0.0
        self.eta_D = 0.0

    def get_LUA_string(self):
        if len(self.systems)==0:
            print "DiatomicRadiator %s has no systems!" % self.name
            sys.exit()
        ostring  = Radiator.get_LUA_string(self)
        ostring += "%s.iTv = %d\n" % ( self.name, self.iTv )
        ostring += "%s.iTr = %d\n" % ( self.name, self.iTr )
        ostring += "%s.I_spin = %d\n" % ( self.name, self.I_spin )
        ostring += "%s.eta_D = %e\n" % ( self.name, self.eta_D )
        ostring += self.level_set.get_LUA_string(mname=self.name)
        ostring += "%s.systems = {\n" % self.name 
        ostring += tab+"systems_list = {"
        for system in self.systems:
            ostring += "'" + system.name + "', "
        ostring += tab+"},\n"
        for system in self.systems:
            ostring += system.get_LUA_string(mname=self.name)
        ostring += "}\n"
        # FIXME: need to check that the level set is compatible with PICS model
        ostring += self.photoionXsection_model.get_LUA_string(species=self.name)
        # FIXME: need to check that the level set is compatible with the QSS model
        if self.E_pop_method=="QSS":
            ostring += self.QSS_model.get_LUA_string(mname=self.name)
        return ostring
        
    def default_data(self):
        self.level_set = self.available_level_sets[self.default_level_set]
        for system in self.available_systems.keys():
            self.systems.append( self.available_systems[system] )
            self.systems[-1].band_set = self.systems[-1].available_band_sets[self.systems[-1].default_band_set]
        self.photoionXsection_model = self.available_photoionXsection_models[self.default_photoionXsection_model]
        self.QSS_model = self.available_QSS_models[self.default_QSS_model]

class DiatomicLevelSet(object):
    """Diatomic level set class"""
    def __init__(self):
        self.levels = []
        self.isp_list = []
        self.comments = "# Description of the diatomic level set"
        
    def get_LUA_string(self, mname):
        ostring  = "%s.level_data = {\n" % ( mname )
        comments = self.comments.replace("#",tab+"--")
        ostring += "%s\n" % comments
        ostring += tab+"n_levels = %d,\n" % ( len(self.levels) )
        ostring += tab+"isp_list = { "
        for isp in self.isp_list:
            ostring += "%d, " % isp
        ostring += "},\n"
        ostring += tab+"-- ===========================================================================================================================================================\n"
        ostring += tab+"--    n      Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  \n"
        ostring += tab+"-- ===========================================================================================================================================================\n"
        for i,level in enumerate(self.levels):
            if i<10: ostring += tab+"ilev_%d  = { %s },\n" % ( i, add_commas(level) )
            else: ostring += tab+"ilev_%d = { %s },\n" % ( i, add_commas(level) )
        ostring += tab+"-- ===========================================================================================================================================================\n"
        ostring += "}\n"
        return ostring
        
class DiatomicSystem:
    """ Diatomic system class"""
    def __init__(self):
        self.name = ""
        self.ie_l = 0
        self.ie_u = 0
        self.band_method = "linebyline"
        self.sigma_nm = 0.5
        self.band_set = DiatomicBandSet()
        self.available_band_sets = {}
        self.comments = "# Description of the diatomic system"
        self.default_band_set = ""
        
    def get_LUA_string(self,mname):
        ostring  = "   %s = {\n" % ( self.name )
        comments = self.comments.replace("#",2*tab+"--")
        ostring += "%s\n" % comments
        ostring += 2*tab+"band_method = '%s',\n" % ( self.band_method )
        ostring += 2*tab+"sigma_nm = %f,\n" % ( self.sigma_nm )
        ostring += 2*tab+"ie_l = %d,\n" % ( self.ie_l )
        ostring += 2*tab+"ie_u = %d,\n" % ( self.ie_u )
        ostring += self.band_set.get_LUA_string(mname=mname,sname=self.name)
        ostring += tab+"},\n"
        return ostring
        
class DiatomicBandSet:
    def __init__(self):
        self.comments = "# Description of the diatomic band set"
        self.uRe_dim = 0
        self.lRe_dim = 0
        self.bands = []
        self.format = "transition moment"
        
    def get_LUA_string(self,mname,sname):
        ostring  = tab*2+"band_data = {\n"
        comments = self.comments.replace("#",tab*3+"--")
        ostring += "%s\n" % comments
        ostring += tab*3+"format = '%s',\n" % ( self.format )
        ostring += tab*3+"uRe_dim = %d,\n" % ( self.uRe_dim )
        ostring += tab*3+"lRe_dim = %d,\n" % ( self.lRe_dim )
        ostring += tab*3+"-- ==========="
        for i in range(self.lRe_dim):
            ostring += "========="
        ostring += "\n"
        ostring += tab*3+"-- Vl =    "
        for i in range(self.lRe_dim):
            ostring += "%d        " % i
        ostring += "\n"
        ostring += tab*3+"-- ==========="
        for i in range(self.lRe_dim):
            ostring += "========="
        ostring += "\n"
        for i,band in enumerate(self.bands):
            if i < 10: ostring += tab*3+"Vu_%d  = { %s },\n" % ( i, add_commas(band) )
            else: ostring += tab*3+"Vu_%d = { %s },\n" % ( i, add_commas(band) )
        ostring += tab*3+"-- ==========="
        for i in range(self.lRe_dim):
            ostring += "========s="
        ostring += "\n"
        ostring += 2*tab+"},\n"
        return ostring

class PolyatomicElectronicLevel:
    def __init__(self, type):
        self.comments = "# Description of the polyatomic electronic level"
        self.type = type
        self.vib_states = []
        self.theta = 0.0
        self.g = 1

    def get_LUA_string(self):
        comments = self.comments.replace('#','       --')
        ostring  = "%s\n" % ( comments )
        ostring += 2*tab + "type = '%s',\n" % self.type
        ostring += 2*tab + "theta = %0.2f,\n" % self.theta
        ostring += 2*tab + "g = %d,\n" % self.g
        ostring += 2*tab + "Nvib_states = %d,\n" % len(self.vib_states)
        ostring += 2*tab + "-- Index               v    g    pv           Gv    Av       Bv       Cv    Dv*10^7  Hv*10^13 Jmax\n"
        for ivs,vs in enumerate(self.vib_states):
            tks = vs[0].split()
            band = tks[0]
            g = int(tks[1])
            Gv = float(tks[2])
            Av = float(tks[3])
            Bv = float(tks[4])
            Cv = float(tks[5])
            Dv = float(tks[6])
            Hv = float(tks[7])
            Jmax = int(tks[8])
            pv = self.calculate_statistical_weight(ivs)
            if 0:
                print "ivs = %d, label = %s, g = %d, pv = %d, Gv = %f" % ( ivs, band, g, pv, Gv )
            if ivs<10:
                spaces = "   "
            elif ivs<100:
                spaces = "  "
            elif ivs<1000:
                spaces = " "
            ostring += 2*tab + "vibstate_%d%s= { '%6s', %2d, %2d, %11.5f, %4.1f, %11.8f, %4.1f, %8.5f, %7.3f, %4d },\n" % ( ivs,spaces,band, g, pv, Gv, Av, Bv, Cv, Dv, Hv, Jmax )
        return ostring
        
    def calculate_statistical_weight(self,ivs):
        vs = self.vib_states[ivs]
        tks = vs[0].split()
        band = tks[0]
        g = int(tks[1])
        if len(band)!=6:
            print "Unexpected size of the vibrational band label!"
            print "Assuming v1=v2=v3=0"
            v = [ 0 ] * 3
        else:
            # currently we only know the format for CO2... (see Lino da Silva phd p218, or Rothman 1981 section 2a)
            v = []
            v.append(int(band[0]))
            v.append(int(band[1]))
            v.append(int(band[3]))
        d = [ 1,2,1 ]
        m = 3
        pv = 1.0
        for im in range(m):
            pv*= factorial(v[im] + d[im] - 1) / factorial(v[im]) / factorial(d[im] - 1)
        return pv
        
class PolyatomicLevelSet(object):
    """Polyatomic level set class"""
    def __init__(self):
        self.levels = []
        self.isp_list = []
        self.comments = "# Description of the polyatomic level set"
        
    def get_LUA_string(self, mname):
        ostring  = "%s.level_data = {\n" % ( mname )
        comments = self.comments.replace("#",tab+"--")
        ostring += "%s\n" % comments
        ostring += tab+"n_levels = %d,\n" % ( len(self.levels) )
        ostring += tab+"isp_list = { "
        for isp in self.isp_list:
            ostring += "%d, " % isp
        ostring += "},\n"
        for i,level in enumerate(self.levels):
            ostring += tab+"ilev_%d  = {\n" % i
            ostring += "%s" % level.get_LUA_string()
            ostring += tab + "},\n"
        ostring += "}\n"
        return ostring

class PolyatomicRadiator(Radiator):
    """Derived diatomic radiator class"""
    def __init__(self, name=""):
        Radiator.__init__(self, name, "polyatomic_radiator")
        self.available_level_sets = {}
        self.photoionXsection_model = PhotoIonXSectionModel()
        self.available_photoionXsection_models = {}
        self.default_photoionXsection_model = ""
        self.iTv = 0
        self.iTr = 0
        self.linear = True
        self.eta_D = 0
        
    def get_LUA_string(self):
        ostring  = Radiator.get_LUA_string(self)
        ostring += "%s.iTv = %d\n" % ( self.name, self.iTv )
        ostring += "%s.iTr = %d\n" % ( self.name, self.iTr )
        ostring += "%s.linear = %d\n" % ( self.name, int(self.linear) )
        ostring += "%s.eta_D = %e\n" % ( self.name, self.eta_D )
        ostring += self.level_set.get_LUA_string(mname=self.name)
        # FIXME: need to check that the level set is compatible with PICS model
        ostring += self.photoionXsection_model.get_LUA_string(species=self.name)
        return ostring

    def default_data(self):
        self.level_set = self.available_level_sets[self.default_level_set]
        self.photoionXsection_model = self.available_photoionXsection_models[self.default_photoionXsection_model]

