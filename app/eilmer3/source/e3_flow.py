"""
e3_flow.py: Functions for reading flow data for 2D and 3D blocks.

This is a grab-bag of classes and functions that help with data handling 
while preparing and post-processing the flow field data.

.. Author: P.A. Jacobs

.. Version: 
     17 March 2008
     15 May 2008 : added reference-function and norm calculations
     01 July 08  : norms returned in a dictionary
     20 Oct 2008 : select surfaces from 3D blocks
"""

from numpy import array, zeros
from gzip import GzipFile
from libprep3 import *
from e3_grid import *
import math
import os
from copy import copy
import struct
import base64
import textwrap

#----------------------------------------------------------------------------
# Utility functions to help set up gas models and flow states in the
# CFCFD codes.

# The FlowCondition class was in the separate file flow_condition2.py
# but it really cannot be used separately from the libgas2 module.
# February 2005: Original code.
# June 2005    : extracted from mbcns_prep.py and eilmer_prep.py into flow_condition2.py
# May 2007     : placed in this file (libgas2.i) as python extension code.
# Jan 2008     : Move the guts to the CFlowCondition class in gas.cxx.
# July 2008    : mu_t, k_t and S added

class FlowCondition(object):
    """
    Python class to organise the setting of each flow condition.

    Now that all of the interesting code has moved to the C++ module gas.cxx,
    this class just provides convenient access from the user's scripts
    as well as some house-keeping needed for writing the config files.
    """
  
    # We will accumulate references to defined objects so that we can
    # later write out an indexed list of flow conditions to the job.config file.
    flowList = []

    __slots__ = 'gmodel', 'flow', 'indx'
    
    def __init__(self, p = 100.0e3, u = 0.0, v = 0.0, w = 0.0,
                 Bx=0.0, By=0.0, Bz=0.0, T = [300.0,],
                 massf = None, label="", tke = 0.0, omega = 1.0,
		 mu_t = 0.0, k_t = 0.0, S = 0, add_to_list=1):
        """
        Create a FlowCondition.

        :param p: static pressure, Pa
        :param u: x-component of velocity, m/s
        :param v: y-component of velocity, m/s
        :param w: z-component of velocity, m/s
        :param Bx: x-component of magnetic field, Tesla
        :param By: y-component of magnetic field, Tesla
        :param Bz: z-component of magnetic field, Tesla
        :param T: list of temperatures (T[0] is static temperature), degrees K
            The length of the list of temperatures must match the number
            of individual energies in the gas model.
        :param massf: mass fractions of the component species
            These may be provided in a number of ways:

            * full list of floats. The length of the list of mass fractions 
                must match the number of species in the previously selected 
                gas model.
            * single float or int that gets used as the first element,
                the rest being set 0.0
            * dictionary of species names with mass fraction values,
                the remainder being set 0.0
            * None provided, results in the first element being 1.0
                and the rest 0.0

	:param label: (optional) string label
        :param tke: turbulence kinetic energy, (m/s)**2
        :param omega: turbulence frequency or pseudo-vorticity, 1/s
	:param mu_t: turbulence viscosity 
	:param k_t: turbulence thermal conductivity
	:param S: shock indicator
 
            * 1 == shock-is-near,
            * 0 == no-shock

        :param add_to_list: flag to indicate that this FlowCondition object 
            should be added to the flowList.  Sometimes we don't want
            to accumulate objects in this list.  For example, when using
            many FlowCondition objects in a user-defined flow evaluation
            function.
        """
	self.gmodel = get_gas_model_ptr()
        if self.gmodel == 0:
            print "ERROR: The gas model has not yet been set."
            print "       A gas model should be set by calling"
            print "       set_type_of_gas(...) before declaring a FlowCondition."
            print "       If unsure, read the documentation."
            sys.exit(-1)
	#
        nsp = self.gmodel.get_number_of_species()
        # species_names = self.gmodel.get_species_names()
        if massf == None:
            massf = [1.0,] + [0.0,]*(nsp-1)
        elif type(massf) is float or type(massf) is int:
            massf = [float(massf),] + [0.0,]*(nsp-1)
        elif type(massf) is dict:
            massf = self.gmodel.to_list(massf)
        massf_sum = 0.0
        for isp in range(nsp): massf_sum += massf[isp]
        if abs(massf_sum - 1.0) > 1.0e-6:
            print "Species mass fractions don't add up to 1.0"
            print "    actual sum=", massf_sum
            print "    massf=", massf
        nmodes = self.gmodel.get_number_of_modes()
        if isinstance(T, list):
            Tlist = T
	elif isinstance(T, float):
            Tlist = [T,] * nmodes
	elif isinstance(T, int):
            Tlist = [float(T),] * nmodes
        else:
            print "Do not know what to do with the supplied value for Tlist."
            print "   T=", T
            print "   type(T)=", type(T)
            raise ValueError, "Bad value supplied for T"
        self.flow = CFlowCondition(self.gmodel, p, u, v, w, Tlist, massf, label, 
				   tke, omega, mu_t, k_t, S, Bx, By, Bz)
        if add_to_list:
            # Indexed FlowCondition objects are later written 
            # to the job.config file.
            self.indx = len(FlowCondition.flowList) # next available index
            FlowCondition.flowList.append(self)
        else:
            # Sometimes, we want to create FlowCondition objects without
            # accumulating them in the flowList.
            self.indx = None
        return

    def __deepcopy__(self, visit):
        """
        Provides a deep copy mechanism for the flow condition.

        It is something like a clone of the original FlowCondition,
        including a containing new CFlowCondition object, but it 
        is not added to the flowList.
        """
        massf = []
        Tlist = []
        nsp = self.gmodel.get_number_of_species()
        nmodes = self.gmodel.get_number_of_modes()
        for isp in range(nsp): massf.append(self.flow.gas.massf[isp])
        for imode in range(nmodes): Tlist.append(self.flow.gas.T[imode])
        newFlow = FlowCondition(p = self.flow.gas.p,
                                u = self.flow.u,
                                v = self.flow.v,
                                w = self.flow.w,
                                Bx = self.flow.Bx,
                                By = self.flow.By,
                                Bz = self.flow.Bz,
                                T = Tlist,
                                massf = massf,
                                label= self.flow.label,
                                tke = self.flow.tke,
                                omega = self.flow.omega,
				mu_t = self.flow.mu_t,
				k_t = self.flow.k_t,
				S = self.flow.S
                                )
        return newFlow

    def __str__(self):
        """
        Produce a string representation that can be used by str() and print.
        """
        str = "FlowCondition("
        str += "p=%g" % (self.flow.gas.p,)
        str += ", u=%g, v=%g, w=%g" % (self.flow.u, self.flow.v, self.flow.w)
        str += ", Bx=%g, By=%g, Bz=%g" % (self.flow.Bx, self.flow.By, self.flow.Bz)
        nsp = self.gmodel.get_number_of_species()
        str += ", massf=["
        for isp in range(nsp): str += "%g," % self.flow.gas.massf[isp]
        str += "]"
        nmodes = self.gmodel.get_number_of_modes()
        str += ", T=["
        for imode in range(nmodes): str += "%g," % self.flow.gas.T[imode]
        str += "]"
        str += ", tke=%g, omega=%g" % (self.flow.tke, self.flow.omega)
        str += ", mu_t=%g, k_t=%g" % (self.flow.mu_t, self.flow.k_t)
        str += ", S=%d" % self.flow.S
        str += ", label=\"" + self.flow.label + "\""
        str += ", add_to_list=%d" % int(self.indx != None)
        str += ")"
        return str

    def write_to_ini_file(self, fp):
        """
        Writes the information to the specified file in .ini format.

        This is used to fill in details in the job.config file.
        """
        fp.write(self.flow.write_to_ini_str(self.indx))
        return

    def to_dict(self):
        """
        Returns the flow data in dictionary form, ready to be written for a cell.
        """
        nsp = self.gmodel.get_number_of_species()
        nmodes = self.gmodel.get_number_of_modes()
        flow_props = {'vel.x':self.flow.u, 'vel.y':self.flow.v, 'vel.z':self.flow.w,
                      'B.x':self.flow.Bx, 'B.y':self.flow.By, 'B.z':self.flow.Bz,
                      'tke':self.flow.tke, 'omega':self.flow.omega,
                      'mu_t':self.flow.mu_t, 'k_t':self.flow.k_t, 'S':self.flow.S}
        flow_props['rho'] = self.flow.gas.rho 
        flow_props['p'] = self.flow.gas.p 
        flow_props['a'] = self.flow.gas.a
        flow_props['mu'] = self.flow.gas.mu
        for imode in range(nmodes):
            flow_props['k[%d]' % imode] = self.flow.gas.k[imode]
        for isp in range(nsp):
            specname = self.gmodel.species_name(isp).replace(' ', '-')
            flow_props['massf[%d]-%s' % (isp,specname)] = self.flow.gas.massf[isp]
        for imode in range(nmodes):
            flow_props['e[%d]' % imode] = self.flow.gas.e[imode]
            flow_props['T[%d]' % imode] = self.flow.gas.T[imode]
        return flow_props

#----------------------------------------------------------------------------
# Reading and writing of Eilmer3-native data files.

def variable_list_for_cell(gdata):
    """
    Returns a list of names for the cell variables that are written
    by the companion function write_cell_data().

    :param gdata: the global-data object (used to control which elements are written)

    .. This function needs to be kept consistent with function write_cell_data(), below,
       and with the corresponding C++ functions in cell.cxx 
       (write_solution_for_cell, read_solution_for_cell and variable_list_for_cell)
    """
    gmodel = get_gas_model_ptr()
    nsp = gmodel.get_number_of_species();
    nmodes = gmodel.get_number_of_modes();
    var_names = ["pos.x", "pos.y", "pos.z", "volume", "rho", "vel.x", "vel.y", "vel.z",]
    if gdata.mhd_flag == 1:
        var_names += ["B.x", "B.y", "B.z"]
    var_names += ["p", "a", "mu"]
    for imode in range(nmodes):
        var_names += [("k[%d]" % imode),]
    var_names += ["mu_t", "k_t", "S"]
    if gdata.radiation_flag == 1: 
        var_names += ["Q_rad_org", "f_rad_org", "Q_rE_rad"]
    var_names.append("tke")
    var_names.append("omega")
    for isp in range(nsp):
        specname = gmodel.species_name(isp).replace(' ', '-')
	var_names.append("massf[%d]-%s" % (isp, specname))
    if nsp > 1: 
        var_names.append("dt_chem")
    for imode in range(nmodes):
	var_names.append("e[%d]" % imode)
	var_names.append("T[%d]" % imode)
    if nmodes > 1: 
        var_names.append("dt_therm")
    return var_names
    
def bgk_list_for_cell(gdata):
    """
    Returns a list of names for the cell bgk variables that are written
    by the companion function write_cell_data().

    :param gdata: the global-data object (used to control which elements are written)
    """
    var_names = ["pos.x", "pos.y", "pos.z", "volume",]
    for gh in range(gdata.velocity_buckets):
        var_names.append("G[%d]" % (gh))
        var_names.append("H[%d]" % (gh))
    return var_names

def quoted_string(vlist):
    """
    Returns a single string with quoted names.
    """
    ost = "\"%s\"" % vlist[0]
    for var in vlist[1:]:
        ost += " \"%s\"" % var
    return ost

def write_cell_data(fp, data, gdata):
    """
    Write the cell data into the specified file (fp).

    :param fp: file object
    :param data: cell data in dictionary form
    :param gdata: the global-data object (used to control which elements are written)

    For Eilmer3 data files, it's all on one line.
    """
    gmodel = get_gas_model_ptr()
    nsp = gmodel.get_number_of_species()
    nmodes = gmodel.get_number_of_modes()
    fp.write("%20.12e %20.12e %20.12e %20.12e" % 
             (data['pos.x'], data['pos.y'], data['pos.z'], data['volume']))
    fp.write(" %20.12e %20.12e %20.12e %20.12e" %
             (data['rho'], data['vel.x'], data['vel.y'], data['vel.z']))
    if gdata.mhd_flag == 1:
        fp.write(" %20.12e %20.12e %20.12e" % (data['B.x'], data['B.y'], data['B.z']))
    fp.write(" %20.12e %20.12e %20.12e" % (data['p'], data['a'], data['mu'],))
    for imode in range(nmodes):
        fp.write(" %20.12e" % data['k[%d]' % imode]) 
    fp.write(" %20.12e %20.12e %1d" % (data['mu_t'], data['k_t'], data['S'],))
    if gdata.radiation_flag == 1:
        fp.write(" %20.12e %20.12e %20.12e" % (0.0,0.0,0.0)) # Zero radiation initially.
    fp.write(" %20.12e %20.12e" % (data['tke'],data['omega']) )
    for isp in range(nsp):
        specname = gmodel.species_name(isp).replace(' ', '-')
        fp.write(" %20.12e" % data['massf[%d]-%s' % (isp, specname)])
    if nsp > 1:
        dt_chem = -1.0
        fp.write(" %20.12e" % dt_chem)
    for imode in range(nmodes):
        fp.write(" %20.12e %20.12e" % (data['e[%d]' % imode], data['T[%d]' % imode],))
    if nmodes > 1:
        dt_therm = -1.0
        fp.write(" %20.12e" % dt_therm)
    fp.write("\n")
    return
    
def write_bgk_data(fp, data, gdata):
    """
    Write the cell data into the specified file (fp).

    :param fp: file object
    :param data: bgk data in dictionary form
    :param gdata: the global-data object (used to control which elements are written)

    For Eilmer3 data files, it's all on one line.
    """
    fp.write("%20.12e %20.12e %20.12e %20.12e" % 
             (data['pos.x'], data['pos.y'], data['pos.z'], data['volume']))
    for gh in range(gdata.velocity_buckets):
        fp.write("%20.12e %20.12e" % (data['G[%d]'%gh], data['H[%d]'%gh]))
    fp.write("\n")
    return

class StructuredGridFlow(object):
    """
    Somewhere to keep the cell and flow data for a structured-grid block.
    """
    def __init__(self):
        self.vars = []
        self.ni = 0
        self.nj = 0
        self.nk = 0
        self.data = {}
        self.gmodel = None
        self.speciesList = []
        self.nsp = 1
        return

    def read(self, fp):
        """
        Read the cell-centre flow data for an entire block, Eilmer3-native format.

        Note that this function cannot cope with spaces inside names.
        """
        buf = fp.readline() # time
        time_stamp = float(buf)
        buf = fp.readline() # variable-name list
        var_list = []
        # We first split the line on spaces ans then strip off quote characters
        # if they are present.  Note that this means that we cannot have spaces
        # inside the names and still expect this simple process to work.
        for token in buf.split():
            self.vars.append(token.strip('"')) # just keep the name
        print self.vars
        buf = fp.readline() # number of cells in each index-direction
        tokens = buf.split()
        ni = int(tokens[0]); self.ni = ni
        nj = int(tokens[1]); self.nj = nj
        nk = int(tokens[2]); self.nk = nk
        self.data = {}
        for variable in self.vars:
            self.data[variable] = zeros((ni,nj,nk),'d')
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    buf = fp.readline() # read cell data
                    tokens = buf.split()
                    for iv in range(len(self.vars)):
                        self.data[self.vars[iv]][i,j,k] = float(tokens[iv])
        return

    def get_cell_data(self, i, j, k, x=0.0, y=0.0, z=0.0, vol=0.0, replace_geom=False):
        """
        Returns the flow data (as a dictionary) for a single cell from a block.
        """
        cell_data = {}
        for variable in self.vars:
            if replace_geom and ( variable == 'pos.x' or variable == 'pos.y' \
                or variable == 'pos.z' or variable == 'volume' ):
                cell_data['pos.x'] = x
                cell_data['pos.y'] = y
                cell_data['pos.z'] = z
                cell_data['volume'] = vol
            else:
                cell_data[variable] = self.data[variable][i,j,k]
        return cell_data
        
    def find_nearest_cell_centre(self, x, y, z):
        """
        Returns the indices of the cell centre nearest the point (x,y,z).
        """
        dx = x - self.data['pos.x'][0,0,0]
        dy = y - self.data['pos.y'][0,0,0]
        dz = z - self.data['pos.z'][0,0,0]
        minDist = math.sqrt(dx*dx + dy*dy + dz*dz)
        imin = 0; jmin = 0; kmin = 0
        for i in range(self.ni):
            for j in range(self.nj):
                for k in range(self.nk):
                    dx = x - self.data['pos.x'][i,j,k]
                    dy = y - self.data['pos.y'][i,j,k]
                    dz = z - self.data['pos.z'][i,j,k]
                    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                    if dist < minDist:
                        minDist = dist
                        imin = i; jmin = j; kmin = k
        return imin, jmin, kmin

    def add_aux_variables(self, cmdLineDict, omegaz=None):
        """
        Adds variables to the data dictionary for each cell in a block.

        We assume a lot about the data that has been read in so,
        we need to skip this function if all is not in place
        """
        for name in ['a', 'rho', 'p', 'vel.x', 'vel.y', 'vel.z', 'e[0]', 'tke']:
            if not name in self.vars:
                print "StructuredGridFlow.add_aux_variables():"
                print "    We are not going to do anything because"
                print "    some essential variables were not present."
                return
        # Always try to attach a gas model, possibly using the command-line argument
        # or by falling back to the expected default file name.
        gmodelFileName = cmdLineDict.get("--gmodel-file", "gas-model.lua")
        print "Attempt to create gas model from file:", gmodelFileName
        if os.path.exists(gmodelFileName):
            self.gmodel = create_gas_model(gmodelFileName)
            self.nsp = self.gmodel.get_number_of_species()
            self.ntm = self.gmodel.get_number_of_modes()
            self.speciesList = [self.gmodel.species_name(isp) for isp in range(self.nsp)]
            print "speciesList=", self.speciesList
        else:
            print "Failed to create gas model."
            self.gmodel = None
            self.nsp = 1
            self.speciesList = []
        # Sift through the command-line dictionary
        # to determine which extra variables should be added.
        add_pitot_p = cmdLineDict.has_key("--add-pitot-p")
        add_total_p = cmdLineDict.has_key("--add-total-p")
        add_total_enthalpy = cmdLineDict.has_key("--add-total-enthalpy")
        add_mach = cmdLineDict.has_key("--add-mach")
        add_molef = cmdLineDict.has_key("--add-molef") and (self.gmodel != None)
        add_trans_coeffs = cmdLineDict.has_key("--add-transport-coeffs") and (self.gmodel != None)
        #
        nic = self.ni; njc = self.nj; nkc = self.nk
        if add_mach:
            self.vars.append("M_local")
            self.data["M_local"] = zeros((nic,njc,nkc),'d')
        if add_pitot_p:
            self.vars.append("pitot_p")
            self.data["pitot_p"] = zeros((nic,njc,nkc),'d')
        if add_total_p:
            self.vars.append("total_p")
            self.data["total_p"] = zeros((nic,njc,nkc),'d')
        if add_total_enthalpy:
            self.vars.append("total_h")
            self.data["total_h"] = zeros((nic,njc,nkc),'d')
        if omegaz is None:
            pass
        else:
            # Add rotating-frame items.
            self.vars.append("M_abs")
            self.data["M_abs"] = zeros((nic,njc,nkc),'d')
            self.vars.append("c.x")
            self.data["c.x"] = zeros((nic,njc,nkc),'d')
            self.vars.append("c.y")
            self.data["c.y"] = zeros((nic,njc,nkc),'d')
            self.vars.append("c.z")
            self.data["c.z"] = zeros((nic,njc,nkc),'d')
            self.vars.append("c.theta")
            self.data["c.theta"] = zeros((nic,njc,nkc),'d')
            self.vars.append("c.r")
            self.data["c.r"] = zeros((nic,njc,nkc),'d')
            self.vars.append("alpha")
            self.data["alpha"] = zeros((nic,njc,nkc),'d')
            self.vars.append("beta")
            self.data["beta"] = zeros((nic,njc,nkc),'d')
        if add_molef:
            for isp in range(self.nsp):
                specname = self.gmodel.species_name(isp).replace(' ', '-')
                varName = "molef[%d]-%s" % (isp, specname)
                self.vars.append(varName)
                self.data[varName] = zeros((nic,njc,nkc),'d')
        if add_trans_coeffs:
            new_vars = [ "sigma", "mu" ]
            for itm in range(self.ntm):
                new_vars.append( "k[%d]"%itm )
            for var in new_vars:
                if var not in self.vars:
                    self.vars.append( var )
                    self.data[var] = zeros((nic,njc,nkc),'d')
        #
        # Now, work through all nodes and add new values.
        for k in range(nkc):
            for j in range(njc):
                for i in range(nic):
                    a = self.data['a'][i,j,k]
                    p = self.data['p'][i,j,k]
                    rho = self.data['rho'][i,j,k]
                    g = a*a*rho/p; # approximation for gamma
                    # Velocity in the block frame of reference.
                    # It may be rotating for turbomachinery calculations.
                    wx = self.data['vel.x'][i,j,k]
                    wy = self.data['vel.y'][i,j,k]
                    wz = self.data['vel.z'][i,j,k]
                    w = math.sqrt(wx*wx + wy*wy + wz*wz)
                    M = w/a
                    if add_mach:
                        self.data['M_local'][i,j,k] = M
                    if add_pitot_p:
                        # Rayleigh Pitot formula
                        if M > 1.0:
                            # Go through the shock and isentropic compression.
                            t1 = (g+1)*M*M/2;
                            t2 = (g+1)/(2*g*M*M - (g-1));
                            pitot_p = p * math.pow(t1,(g/(g-1))) * math.pow(t2,(1/(g-1)));
                        else:
                            # Isentropic compression only.
                            t1 = 1 + 0.5*(g-1)*M*M;
                            pitot_p = p * math.pow(t1,(g/(g-1)));
                        self.data['pitot_p'][i,j,k] = pitot_p
                    if add_total_p:
                        # Isentropic process only.
                        t1 = 1 + 0.5*(g-1)*M*M;
                        total_p = p * math.pow(t1,(g/(g-1)));
                        self.data['total_p'][i,j,k] = total_p
                    if add_total_enthalpy:
                        e0 = self.data['e[0]'][i,j,k]
                        tke = self.data['tke'][i,j,k]
                        # Sum up the bits of energy,
                        # forgetting the multiple energy modes, for the moment.
                        total_h = p/rho + e0 + 0.5*w*w + tke;
                        self.data['total_h'][i,j,k] = total_h
                    if omegaz is None:
                        pass
                    else:
                        # Add the rotating-frame variables.
                        x = self.data['pos.x'][i,j,k]
                        y = self.data['pos.y'][i,j,k]
                        r = math.sqrt(x*x + y*y)
                        sintheta = y/r
                        costheta = x/r
                        # The following should get us back to the lab (or absolute) frame.
                        cx = wx - omegaz*y
                        cy = wy + omegaz*x
                        cz = wz
                        ctheta = -sintheta*cx + costheta*cy
                        cr = costheta*cx + sintheta*cy
                        c = math.sqrt(cr*cr + ctheta*ctheta + cz*cz)
                        alpha = math.atan2(ctheta,cz)
                        beta = math.asin(cr/c)
                        M_abs = c/a
                        self.data['M_abs'][i,j,k] = M_abs
                        self.data['c.x'][i,j,k] = cx
                        self.data['c.y'][i,j,k] = cy
                        self.data['c.z'][i,j,k] = cz
                        self.data['c.theta'][i,j,k] = ctheta
                        self.data['c.r'][i,j,k] = cr
                        self.data['alpha'][i,j,k] = alpha
                        self.data['beta'][i,j,k] = beta
                    if add_molef:
                        # Accumulate the mass fractions and then compute the mole fractions
                        # in the context of the gas model.
                        massf = []
                        for isp in range(self.nsp):
                            specname = self.gmodel.species_name(isp).replace(' ', '-')
                            massf.append(self.data['massf[%d]-%s' % (isp, specname)][i,j,k])
                        # molef = convert_massf2molef(massf, self.gmodel.M()) # works OK
                        molef = self.gmodel.to_molef(massf)
                        for isp in range(self.nsp):
                            specname = self.gmodel.species_name(isp).replace(' ', '-')
                            self.data['molef[%d]-%s' % (isp, specname)][i,j,k] = molef[isp]
                    if add_trans_coeffs:
                        # recompute transport properties with the provided gas-model file
                        Q = Gas_data(self.gmodel)
                        nsp = self.gmodel.get_number_of_species()
                        nmodes = self.gmodel.get_number_of_modes()
                        Q.rho = self.data['rho'][i,j,k]
                        for isp in range(nsp):
                            specname = self.gmodel.species_name(isp).replace(' ', '-')
                            Q.massf[isp] = self.data['massf[%d]-%s' % (isp, specname)][i,j,k]
                        for imode in range(nmodes):
                            Q.T[imode] = self.data['T[%d]' % imode][i,j,k]
                        self.gmodel.eval_thermo_state_rhoT(Q)
                        self.gmodel.eval_transport_coefficients(Q)
                        self.data['sigma'][i,j,k] = Q.sigma
                        self.data['mu'][i,j,k] = Q.mu
                        for imode in range(nmodes):
                            self.data['k[%d]'%imode][i,j,k] = Q.k[imode]
                            
        # end of adding new data values for a block
        return

    def write_gnuplot_header(self, fp):
        """
        Write a comment line that is compatible with GNUPlot.
        """
        fp.write("# Variables:")
        for i, var in enumerate(self.vars):
            # i+1 because Python numbers from 0, but GNUPlot
            # starts numbering columns from 1
            fp.write(" %d:%s" % (i+1, var))
        fp.write("\n")
        return

    def write_gnuplot_data_for_cell(self, fp, i, j, k):
        """
        Write the flow data for a single cell to the specified file.
        """
        for var in self.vars:
            # print "var=", var
            fp.write(" %e" % self.data[var][i,j,k])
        fp.write("\n")
        return

# end of class StructuredGridFlow

def read_time_from_flow_file(rootName, tindx, zipFiles=False):
    """
    We'll find the simulation time on the first lins of the flow file.
    """
    # tindx may be an integer, or already a string such as "xxxx"
    if type(tindx) is int:
        tindx_str = "%04d" % tindx
    elif type(tindx) is str:
        tindx_str = tindx
    else:
        raise RuntimeException("WTF: tindx is neither an int nor string.")
    fileName = rootName+".flow"+(".b0000.t%04s" % tindx_str)
    fileName = os.path.join("flow", "t%04s" % tindx_str, fileName)
    print "Read simulation time from", fileName
    if zipFiles: 
        fp = GzipFile(fileName+".gz", "rb")
    else:
        fp = open(fileName, "r")
    line = fp.readline().strip().strip("\0")
    t = float(line)
    fp.close()
    return t

def read_all_blocks(rootName, nblock, tindx, zipFiles=False, movingGrid=False):
    """
    Returns all grids and flow blocks for a single flow solution.

    """
    # tindx may be an integer, or already a string such as "xxxx"
    if type(tindx) is int:
        tindx_str = "%04d" % tindx
    elif type(tindx) is str:
        tindx_str = tindx
    else:
        raise RuntimeException("WTF: tindx is neither an int nor string.")
    grid = []; flow = []
    for jb in range(nblock):
        if not movingGrid:
            fileName = rootName+".grid"+(".b%04d.t%04s" % (jb, "0000"))
            fileName = os.path.join("grid", "t0000", fileName)
        else:
            fileName = rootName+".grid"+(".b%04d.t%04s" % (jb, tindx_str))
            fileName = os.path.join("grid", "t%04s" % tindx_str, fileName)
        print "Read cell-vertex data from", fileName
        grid.append(StructuredGrid())
        if zipFiles: 
            fp = GzipFile(fileName+".gz", "rb")
        else:
            fp = open(fileName, "r")
        grid[-1].read(fp)
        fp.close()
        #
        fileName = rootName+".flow"+(".b%04d.t%04s" % (jb, tindx_str))
        fileName = os.path.join("flow", "t%04s" % tindx_str, fileName)
        print "Read cell-centre flow data from", fileName
        if zipFiles: 
            fp = GzipFile(fileName+".gz", "rb")
        else:
            fp = open(fileName, "r")
        flow.append(StructuredGridFlow())
        flow[-1].read(fp)
        fp.close()
        #
        fileName = rootName+".BGK"+(".b%04d.t%04s" % (jb, tindx_str))
        fileName = os.path.join("flow", "t%04s" % tindx_str, fileName)
        # only open a BGK file if one is there
        if os.path.exists(fileName) | os.path.exists(fileName+".gz"):
            print "Read velocity distribution data from", fileName
            if zipFiles: 
                fp = GzipFile(fileName+".gz", "rb")
            else:
                fp = open(fileName, "r")
            bgk = StructuredGridFlow()
            bgk.read(fp)
            fp.close()
            print "Append velocity distribution data to flow data"
            for var in bgk.vars:
                if var not in flow[-1].vars:
                    flow[-1].vars.append(var)
                    flow[-1].data[var] = bgk.data[var]
    #
    if grid[0].nk == 1:
        dimensions = 2
    else:
        dimensions = 3
    return grid, flow, dimensions

def add_auxiliary_variables(nblock, flow, cmdLineDict, omegaz_list=None):
    """
    Adds variables to the data dictionary for each cell in each block.
    """
    for jb in range(nblock):
        if omegaz_list == None:
            omegaz = None
        else:
            omegaz = omegaz_list[jb]
        flow[jb].add_aux_variables(cmdLineDict, omegaz)
    return

def locate_cell_and_block(grid, flow, dimensions, 
                          i_found, j_found, k_found, jb_found,
                          x, y, z=0.0):
    """
    Returns the indices that select a particular cell that contains 
    point (x,y,z) or has cell-centre nearest the point, 
    together with the containing block of flow_data.

    Can be used for interpolating flow data from within another solution.
    
    .. 10.03.2010 - Modified to use suggest_better_cell() function, DFP
    .. 25.04.2010 - moved suggest_better_cell() to class StructuredGrid, PJ.
    """
    # Step 1 : Starts by checking if point is in the last found block
    jb = jb_found
    di = 1; dj = 1; dk = 1
    i = i_found; j = j_found; k = k_found
    while 1:
    	di, dj, dk = grid[jb].suggest_better_cell(i, j, k, x, y, z, dimensions)
    	if abs(di) + abs(dj) + abs(dk) == 0:
    	    return jb, i, j, k
    	i += di; j += dj; k += dk
    	if i >= flow[jb].ni or i<0: break
    	if j >= flow[jb].nj or j<0: break
    	if k >= flow[jb].nk or k<0: break
    # Step 2: Try the other blocks ahead of jb_found
    for jb in range(jb_found+1, len(flow)):
    	di = 1; dj = 1; dk = 1
    	i = 0; j = 0; k = 0
    	while 1:
    	    di, dj, dk = grid[jb].suggest_better_cell(i, j, k, x, y, z, dimensions)
    	    if abs(di) + abs(dj) + abs(dk) == 0:
    	    	return jb, i, j, k
    	    i += di; j += dj; k += dk
    	    if i >= flow[jb].ni or i<0: break
    	    if j >= flow[jb].nj or j<0: break
    	    if k >= flow[jb].nk or k<0: break
    # Step 3: Try the other blocks before jb_found
    for jb in range(0, jb_found):
    	di = 1; dj = 1; dk = 1
    	i = 0; j = 0; k = 0
    	while 1:
    	    di, dj, dk = grid[jb].suggest_better_cell(i, j, k, x, y, z, dimensions)
    	    if abs(di) + abs(dj) + abs(dk) == 0:
    	    	return jb, i, j, k
    	    i += di; j += dj; k += dk
    	    if i >= flow[jb].ni or i<0: break
    	    if j >= flow[jb].nj or j<0: break
    	    if k >= flow[jb].nk or k<0: break
    # Step 4: If the search is still unsuccessful at this point,
    #         we should try to find a nearest neighbouring cell.
    #         Or maybe the geometries of both grids are not similar.
    px=flow[jb_found].data['pos.x'][i_found,j_found,k_found]
    py=flow[jb_found].data['pos.y'][i_found,j_found,k_found]
    pz=flow[jb_found].data['pos.z'][i_found,j_found,k_found]
    dist = math.sqrt((x-px)**2 + (y-py)**2 + (z-pz)**2)
    best_block_indx = jb_found; best_i = i_found; best_j = j_found; best_k = k_found; min_dist = dist
    for jb in range(0, len(flow)):
        for k in range(flow[jb].nk):
            for j in range(flow[jb].nj):
                for i in range(flow[jb].ni):
                    px=flow[jb].data['pos.x'][i,j,k]
                    py=flow[jb].data['pos.y'][i,j,k]
                    pz=flow[jb].data['pos.z'][i,j,k]
                    dist = math.sqrt((x-px)**2 + (y-py)**2 + (z-pz)**2)
                    if dist < min_dist:
                        min_dist = dist; best_block_indx = jb
                        best_i = i; best_j = j; best_k = k
    # If we have reached this point, then we have not succeeded 
    # in locating the cell for the specified point, so we return
    # with the closest cell.
    return best_block_indx, best_i, best_j, best_k


class ExistingSolution(object):
    """
    A place to store an existing solution and to provide interpolated flow data.

    A common use case is to transfer a previously computed solution into
    the initial conditions of a new simulation.
    
    The interpolate_flow_condition() method is the primary way to access the 
    flow data when e3prep.py is filling in the initial flow state for a block.
    This function returns the flow data in a form suitable for write_cell_data()
    and can be used in the user's setup script for e3prep.py.
    """
    def __init__(self, rootName, solutionWorkDir, nblock, tindx, 
                 dimensions=2, assume_same_grid=0, zipFiles=1,
                 add_velocity=Vector(0.0,0.0,0.0)):
        """
        Reads and stores an existing solution.

        :param rootName: job name that will be used to build file names
        :param solutionWorkDir: the directory where we'll find our ExistingSolution files.
        :param nblock: number of blocks in the ExistingSolution data set
        :param tindx: the time index to select 0..9999:
            Do not specify with leading zeros because the Python interpreter
            will assume that you want to count the time index in octal.

        :param dimensions: number of spatial dimensions for the ExistingSolution
        :param assume_same_grid: decide how to locate corresponding cells
 
            * 0 == searches for corresponding cells
              As Rainer found, this can be agonisingly slow for large grids.
            * 1 == omits the search for the corresponding cell
              For the impatient.

        :param zipFiles: to use gzipped files, or not
        :param add_velocity: is the value to be aded to the old-solution's velocity.
            If it is a float value, it becomes the x-component. Otherwise, 
            a Vector value should be supplied and all components will be used.
        """
        currentWorkDir = os.getcwd()  # Remembers the current working directory
        os.chdir(solutionWorkDir)     # Changes to the working directory where
                                      # the existing solution is stored
        self.rootName = rootName # remember this name for later reporting
        if type(add_velocity) is float:
            # Assume that a single component relates to the x-direction.
            self.add_velocity = Vector(add_velocity, 0.0, 0.0)
        else:
            try:
                self.add_velocity = Vector(add_velocity.x, add_velocity.y, add_velocity.z)
            except:
                print "ExistingSolution expected a float or Vector for add_velocity."
                print "    received a type", type(add_velocity), "with value", add_velocity
                raise TypeError("add_velocity parameter is incorrect")
        self.dimensions = dimensions
        self.assume_same_grid = assume_same_grid
        self.grid, self.flow, stored_dims = read_all_blocks(rootName, nblock, tindx, zipFiles)
        if stored_dims != dimensions:
            print "Oops: store_dims=", stored_dims, "dimensions=", dimensions
            print "These values of dimensions should be the same."
        os.chdir(currentWorkDir)   # Changes back to the original working directory
        return

    def interpolate_flow_condition(self, x, y, z, vol, i, j, k, blk_indx):
        """
        Returns the flow data for a given location.
        """
        if self.assume_same_grid:
            # no need to search because we assume i,j,k correspond to correct location
            block_flow = self.flow[blk_indx]
            block_number = blk_indx
            if self.bgk:
                block_bgk = self.bgk[blk_indx]
        else:
            # we need to locate the appropriate cell
            jb , i, j, k = locate_cell_and_block(self.grid, self.flow, self.dimensions,
                                                 i, j, k, blk_indx, x, y, z)
            block_flow = self.flow[jb]
            block_number = jb
            if False and self.bgk:
                block_bgk = self.bgk[jb]
                
        flow_data_for_cell = block_flow.get_cell_data(i, j, k, x, y, z, vol, replace_geom=True)
        flow_data_for_cell['vel.x'] += self.add_velocity.x
        flow_data_for_cell['vel.y'] += self.add_velocity.y
        flow_data_for_cell['vel.z'] += self.add_velocity.z
        
        if False and self.bgk:
            bgk_data_for_cell = block_bgk.get_cell_data(i, j, k, x, y, z, vol, replace_geom=True)
            return (flow_data_for_cell, bgk_data_for_cell, i, j, k, block_number)
        else:
            return (flow_data_for_cell, i, j, k, block_number)

#---------------------------------------------------------------------------------
# VTK-related functions

def write_VTK_XML_unstructured_file(fp, grid, flow, binary_format):
    """
    Write the cell-centred flow data from a single block 
    as an unstructured grid of finite-volume cells.

    :param fp: reference to a file object
    :param grid: single-block grid of vertices
    :param flow: single-block of cell-centre flow data
    :param binary_format: if True, most of the data will be written as raw binary 
        (appended) else it will be written as ascii (in place).
        We have not used base64 encoding; too much pain with Paraview.
    """
    if binary_format:
        binary_data_string = ""
        binary_data_offset = 0
    niv = grid.ni; njv = grid.nj; nkv = grid.nk
    nic = flow.ni; njc = flow.nj; nkc = flow.nk
    two_D = (nkv == 1)
    NumberOfPoints = niv * njv * nkv
    NumberOfCells = nic * njc * nkc
    fp.write("<VTKFile type=\"UnstructuredGrid\" byte_order=\"BigEndian\">\n")
    fp.write("<UnstructuredGrid>")
    fp.write("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n" %
             (NumberOfPoints, NumberOfCells))
    #
    fp.write("<Points>\n")
    fp.write(" <DataArray type=\"Float32\" NumberOfComponents=\"3\"")
    if binary_format:
        fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
        binary_data = r""
    else:
        fp.write(" format=\"ascii\">\n")
    vtx_number = 0
    vtx_id = {}
    for k in range(nkv):
        for j in range(njv):
            for i in range(niv):
                vtx_id[(i,j,k)] = vtx_number
                x,y,z = uflowz(grid.x[i,j,k]), uflowz(grid.y[i,j,k]), uflowz(grid.z[i,j,k])
                if binary_format:
                    binary_data += struct.pack('> f f f', x, y, z)
                else:
                    fp.write(" %e %e %e\n" % (x,y,z))
                vtx_number += 1
    fp.write(" </DataArray>\n")
    if binary_format:
        binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
        binary_data_string += binary_data_count
        binary_data_string += binary_data
        binary_data_offset += len(binary_data_count) + len(binary_data)
    fp.write("</Points>\n")
    #
    fp.write("<Cells>\n")
    fp.write(" <DataArray type=\"Int32\" Name=\"connectivity\"")
    if binary_format:
        fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
        binary_data = r""
    else:
        fp.write(" format=\"ascii\">\n")
    for k in range(nkc):
        for j in range(njc):
            for i in range(nic):
                if two_D:
                    ids = tuple([vtx_id[(i,j,k)], vtx_id[(i+1,j,k)],
                                 vtx_id[(i+1,j+1,k)], vtx_id[(i,j+1,k)]])
                    if binary_format:
                        binary_data += struct.pack('> i i i i', ids[0], ids[1], ids[2], ids[3])
                    else:
                        fp.write(" %d %d %d %d\n" % ids)
                else:
                    ids = tuple([vtx_id[(i,j,k)], vtx_id[(i+1,j,k)], 
                                 vtx_id[(i+1,j+1,k)], vtx_id[(i,j+1,k)],
                                 vtx_id[(i,j,k+1)], vtx_id[(i+1,j,k+1)], 
                                 vtx_id[(i+1,j+1,k+1)], vtx_id[(i,j+1,k+1)]])
                    if binary_format:
                        binary_data += struct.pack('> i i i i i i i i', ids[0], ids[1], ids[2],
                                                   ids[3], ids[4], ids[5], ids[6], ids[7])
                    else:
                        fp.write(" %d %d %d %d %d %d %d %d\n" % ids)
    fp.write(" </DataArray>\n")
    if binary_format:
        binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
        binary_data_string += binary_data_count
        binary_data_string += binary_data
        binary_data_offset += len(binary_data_count) + len(binary_data)
    #
    fp.write(" <DataArray type=\"Int32\" Name=\"offsets\"")
    if binary_format:
        fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
        binary_data = r""
    else:
        fp.write(" format=\"ascii\">\n")
    # Since all of the point-lists are concatenated, these offsets into the connectivity
    # array specify the end of each cell.
    for k in range(nkc):
        for j in range(njc):
            for i in range(nic):
                if two_D:
                    conn_offset = 4*(1+i+j*nic)
                else:
                    conn_offset = 8*(1+i+j*nic+k*(nic*njc))
                if binary_format:
                    binary_data += struct.pack('> i', conn_offset)
                else:
                    fp.write(" %d\n" % conn_offset)
    fp.write(" </DataArray>\n")
    if binary_format:
        binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
        binary_data_string += binary_data_count
        binary_data_string += binary_data
        binary_data_offset += len(binary_data_count) + len(binary_data)
    #
    fp.write(" <DataArray type=\"UInt8\" Name=\"types\"")
    if binary_format:
        fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
        binary_data = r""
    else:
        fp.write(" format=\"ascii\">\n")
    if two_D:
        type_value = 9 # VTK_QUAD
    else:
        type_value = 12 # VTK_HEXAHEDRON
    for k in range(nkc):
        for j in range(njc):
            for i in range(nic):
                if binary_format:
                    binary_data += struct.pack('> B', type_value)
                else:
                    fp.write(" %d\n" % type_value)
    fp.write(" </DataArray>\n")
    if binary_format:
        binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
        binary_data_string += binary_data_count
        binary_data_string += binary_data
        binary_data_offset += len(binary_data_count) + len(binary_data)
    fp.write("</Cells>\n")
    #
    fp.write("<CellData>\n")
    # Write variables from the dictionary.
    for var in flow.vars:
        fp.write(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"" % (var))
        if binary_format:
            fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
            binary_data = r""
        else:
            fp.write(" format=\"ascii\">\n")
        for k in range(nkc):
            for j in range(njc):
                for i in range(nic):
                    if binary_format:
                        binary_data += struct.pack('> f', uflowz(flow.data[var][i,j,k]))
                    else:
                        fp.write(" %e\n" % uflowz(flow.data[var][i,j,k]))
        fp.write(" </DataArray>\n")
        if binary_format:
            binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
            binary_data_string += binary_data_count
            binary_data_string += binary_data
            binary_data_offset += len(binary_data_count) + len(binary_data)
    #
    # Write the special variables:
    # i.e. variables constructed from those in the dictionary.
    fp.write(" <DataArray Name=\"vel.vector\" type=\"Float32\" NumberOfComponents=\"3\"")
    if binary_format:
        fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
        binary_data = r""
    else:
        fp.write(" format=\"ascii\">\n")
    for k in range(nkc):
        for j in range(njc):
            for i in range(nic):
                x, y, z = (uflowz(flow.data['vel.x'][i,j,k]),
                           uflowz(flow.data['vel.y'][i,j,k]),
                           uflowz(flow.data['vel.z'][i,j,k]))
                if binary_format:
                    binary_data += struct.pack('> f f f', x, y, z)
                else:
                    fp.write(" %e %e %e\n" % (x,y,z))
    fp.write(" </DataArray>\n")
    if binary_format:
        binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
        binary_data_string += binary_data_count
        binary_data_string += binary_data
        binary_data_offset += len(binary_data_count) + len(binary_data)
    #
    if 'c.x' in flow.vars:
        fp.write(" <DataArray Name=\"c.vector\" type=\"Float32\" NumberOfComponents=\"3\"")
        if binary_format:
            fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
            binary_data = r""
        else:
            fp.write(" format=\"ascii\">\n")
        for k in range(nkc):
            for j in range(njc):
                for i in range(nic):
                    x, y, z = (uflowz(flow.data['c.x'][i,j,k]), 
                               uflowz(flow.data['c.y'][i,j,k]), 
                               uflowz(flow.data['c.z'][i,j,k]))
                    if binary_format:
                        binary_data += struct.pack('> f f f', x, y, z)
                    else:
                        fp.write(" %e %e %e\n" % (x, y, z) )
        fp.write(" </DataArray>\n")
        if binary_format:
            binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
            binary_data_string += binary_data_count
            binary_data_string += binary_data
            binary_data_offset += len(binary_data_count) + len(binary_data)
    #
    if 'B.x' in flow.vars:
        fp.write(" <DataArray Name=\"B.vector\" type=\"Float32\" NumberOfComponents=\"3\">\n")
        for k in range(nkc):
            for j in range(njc):
                for i in range(nic):
                    x, y, z = (uflowz(flow.data['B.x'][i,j,k]), 
                               uflowz(flow.data['B.y'][i,j,k]), 
                               uflowz(flow.data['B.z'][i,j,k]))
                    if binary_format:
                        binary_data += struct.pack('> f f f', x, y, z)
                    else:
                        fp.write(" %e %e %e\n" % (x, y, z) )
        fp.write(" </DataArray>\n")
        if binary_format:
            binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
            binary_data_string += binary_data_count
            binary_data_string += binary_data
            binary_data_offset += len(binary_data_count) + len(binary_data)
    #
    fp.write("</CellData>\n")
    fp.write("</Piece>\n")
    fp.write("</UnstructuredGrid>\n")
    if binary_format:
        fp.write("<AppendedData encoding=\"raw\">\n")
        fp.write('_'+binary_data_string)
        fp.write("</AppendedData>\n")
    fp.write("</VTKFile>\n")
    return

def begin_Visit_file(rootName, nblock):
    # Will be handy to have a Visit file, also.
    # For each time index, this justs lists the names of the files for individual blocks.
    plotPath = "plot"
    if not os.access(plotPath, os.F_OK):
        os.makedirs(plotPath)
    fileName = rootName+".visit"
    fileName = os.path.join(plotPath, fileName)
    visitFile = open(fileName, "w")
    visitFile.write("!NBLOCKS %d\n" % nblock)
    visitFile.close()
    return

def begin_PVD_file(rootName):
    # Will be handy to have a Paraview collection file, also.
    # For each time index, this justs lists the name of the top-level .pvtu file.
    plotPath = "plot"
    if not os.access(plotPath, os.F_OK):
        os.makedirs(plotPath)
    fileName = rootName+".pvd"
    fileName = os.path.join(plotPath, fileName)
    pvdFile = open(fileName, "w")
    pvdFile.write("<?xml version=\"1.0\"?>\n")
    pvdFile.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
    pvdFile.write("<Collection>\n")
    pvdFile.close()
    return

def finish_PVD_file(rootName):
    plotPath = "plot"
    if not os.access(plotPath, os.F_OK):
        os.makedirs(plotPath)
    fileName = rootName+".pvd"
    fileName = os.path.join(plotPath, fileName)
    pvdFile = open(fileName, "a")
    pvdFile.write("</Collection>\n")
    pvdFile.write("</VTKFile>\n")
    pvdFile.close()
    return

def write_VTK_XML_files(rootName, tindx, nblock, grid, flow, t, binary_format):
    """
    Writes the top-level (coordinating) parallel/partitioned-VTK file.

    :param rootName: specific file names are built by adding bits to this name
    :param tindx: integer or string such as "xxxx"
    :param nblock: integer
    :param grid: list of StructuredGrid objects
    :param flow: list of StructuredGridFlow objects
    :param t: simulation time corresponding to the flow data at this tindx 
    :param binary_format: if True, most of the data will be written as raw binary 
        (appended) else it will be written as ascii (in place).
    """
    plotPath = "plot"
    if not os.access(plotPath, os.F_OK):
        os.makedirs(plotPath)
    # tindx may be an integer, or already a string such as "xxxx"
    if type(tindx) is int:
        tindx_str = "%04d" % tindx
    elif type(tindx) is str:
        tindx_str = tindx
    else:
        raise RuntimeException("WTF: tindx is neither an int nor string.")
    fileName = rootName+(".t%04s" % tindx_str)+".pvtu"
    fileName = os.path.join(plotPath, fileName)
    pvtuFile = open(fileName, "w")
    pvtuFile.write("<VTKFile type=\"PUnstructuredGrid\">\n")
    pvtuFile.write("<PUnstructuredGrid GhostLevel=\"0\">")
    pvtuFile.write("<PPoints>\n")
    pvtuFile.write(" <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n")
    pvtuFile.write("</PPoints>\n")
    pvtuFile.write("<PCellData>\n")
    for var in flow[0].vars:
        pvtuFile.write(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"/>\n" % (var))
    pvtuFile.write(" <PDataArray Name=\"vel.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n")
    if 'c.x' in flow[0].vars:
        pvtuFile.write(" <PDataArray Name=\"c.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n")
    if 'B.x' in flow[0].vars:
        pvtuFile.write(" <PDataArray Name=\"B.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n")
    pvtuFile.write("</PCellData>\n")
    for jb in range(nblock):
        fileName = rootName+(".b%04d.t%04s" % (jb, tindx_str))+".vtu"
        # We write the short version of the fileName into the pvtu file.
        pvtuFile.write("<Piece Source=\"%s\"/>\n" % fileName)
        # but use the long version to actually open it.
        fileName = os.path.join(plotPath, fileName)
        vtuFile = open(fileName, "wb") # We may be writing some binary data.
        write_VTK_XML_unstructured_file(vtuFile, grid[jb], flow[jb], binary_format)
        vtuFile.close()
    pvtuFile.write("</PUnstructuredGrid>\n")
    pvtuFile.write("</VTKFile>\n")
    pvtuFile.close()
    # Will be handy to have a Visit file, also.
    # This justs lists the names of the files for individual blocks.
    fileName = rootName+".visit"
    fileName = os.path.join(plotPath, fileName)
    # Note that we append to the visit file for each tindx.
    visitFile = open(fileName, "a")
    for jb in range(nblock):
        fileName = rootName+(".b%04d.t%04s" % (jb, tindx_str))+".vtu"
        visitFile.write("%s\n" % fileName)
    visitFile.close()
    # Will be handy to have a Paraview PVD file, also.
    # This justs lists the top-level .pvtu files.
    fileName = rootName+".pvd"
    fileName = os.path.join(plotPath, fileName)
    # Note that we append to the .pvd file for each tindx.
    pvdFile = open(fileName, "a")
    fileName = rootName+(".t%04s" % (tindx_str,))+".pvtu"
    pvdFile.write("<DataSet timestep=\"%e\" group=\"\" part=\"0\" file=\"%s\"/>\n" % (t, fileName))
    pvdFile.close()
    return

#---------------------------------------------------------------------------------
# Tecplot-related functions

def write_Tecplot_file(rootName, tindx, nblock, grid, flow, t):
    """
    Write the TECPLOT ASCII file in BLOCK form, CELLCENTERED flow data.

    :param rootName: specific file names are built by adding bits to this name
    :param tindx: integer
    :param nblock: integer
    :param grid: list of StructuredGrid objects
    :param flow: list of StructuredGridFlow objects
    :param t: float value for the time that gets written into the file

    Note that the Tecplot360 ASCII data format has a limitation of 32000 characters 
    on a line.  This may be a problem for the simple code below if we exceed nni=1600.
    For the moment, retain the simpler code rather than making it cope with nni>1600.
    """
    plotPath = "plot"
    if not os.access(plotPath, os.F_OK):
        os.makedirs(plotPath)
    # tindx may be an integer, or already a string such as "xxxx"
    if type(tindx) is int:
        tindx_str = "%04d" % tindx
    elif type(tindx) is str:
        tindx_str = tindx
    else:
        raise RuntimeException("WTF: tindx is neither an int nor string.")
    fileName = rootName+(".t%04s" % tindx_str)+".tec"
    fileName = os.path.join(plotPath, fileName)
    fp = open(fileName, "w")
    fp.write("TITLE=\"Job=%s time=%e\"\n" % (rootName, t))
    fp.write("VARIABLES= \"X\", \"Y\", \"Z\"")
    n_centered_vars = 0
    for var in flow[0].vars:
        if var == "pos.x" or var == "pos.y" or var == "pos.z": continue
        fp.write(", \"%s\"" % (var))
        n_centered_vars += 1
    fp.write("\n")
    centered_list_str = str(range(4,4+n_centered_vars))
    for jb in range(nblock):
        nic = flow[jb].ni; njc = flow[jb].nj; nkc = flow[jb].nk
        niv = grid[jb].ni; njv = grid[jb].nj; nkv = grid[jb].nk
        fp.write("ZONE I=%d J=%d K=%d DATAPACKING=BLOCK" % (niv,njv,nkv,))
        fp.write(" SOLUTIONTIME=%e" % t)
        fp.write(" VARLOCATION=(%s=CELLCENTERED) T=\"block-%d\"\n" % (centered_list_str,jb,))
        fp.write("# cell-vertex pos.x\n")
        for k in range(nkv):
            for j in range(njv):
                for i in range(niv):
                    fp.write(" %e" % uflowz(grid[jb].x[i,j,k]))
                fp.write("\n")
        fp.write("# cell-vertex pos.y\n")
        for k in range(nkv):
            for j in range(njv):
                for i in range(niv):
                    fp.write(" %e" % uflowz(grid[jb].y[i,j,k]))
                fp.write("\n")
        fp.write("# cell-vertex pos.z\n")
        for k in range(nkv):
            for j in range(njv):
                for i in range(niv):
                    fp.write(" %e" % uflowz(grid[jb].z[i,j,k]))
                fp.write("\n")
        for var in flow[jb].vars:
            if var == "pos.x" or var == "pos.y" or var == "pos.z": continue
            fp.write("# cell-centre %s\n" % var)
            for k in range(nkc):
                for j in range(njc):
                    for i in range(nic):
                        fp.write(" %e" % uflowz(flow[jb].data[var][i,j,k]))
                    fp.write("\n")
    fp.close()
    return

#-----------------------------------------------------------------------------
# Plot3D related functions
def write_plot3d_files(rootName, tindx, nblock, grid, flow, t):
    """
    Write the Plote3D files.

    :param rootName: specific file names are built by adding bits to this name
    :param tindx: integer
    :param nblock: integer
    :param grid: list of StructuredGrid objects
    :param flow: list of StructuredGridFlow objects
    :param t: float value for time (that gets ignored)
    """
    plotPath = "plot"
    if not os.access(plotPath, os.F_OK):
        os.makedirs(plotPath)
    # tindx may be an integer, or already a string such as "xxxx"
    if type(tindx) is int:
        tindx_str = "%04d" % tindx
    elif type(tindx) is str:
        tindx_str = tindx
    else:
        raise RuntimeException("WTF: tindx is neither an int nor string.")
    #
    # 1. write the 'grid' file.
    #    This is just cell-centre locations
    gfileName = rootName+(".t%04s" % tindx_str)+".g"
    gfileName = os.path.join(plotPath, gfileName)
    fp = open(gfileName, 'w')
    fp.write("%d\n" % len(flow))
    for f in flow:
        ni = f.ni; nj = f.nj; nk = f.nk
        if nk == 1:
            fp.write("%d %d\n" % (ni, nj))
        else:
            fp.write("%d %d %d\n" % (ni, nj, nk))

        coords = ["pos.x", "pos.y", "pos.z"]
        #
        for c in coords:
            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        fp.write("%20.12e\n" % f.data[c][i,j,k])
    #
    # 2. Write the .nam file
    #    We remove pos.x, pos.y, pos.z from the list of variables
    nfileName = rootName+(".t%04s" % tindx_str)+".nam"
    nfileName = os.path.join(plotPath, nfileName)
    fp = open(nfileName, 'w')
    for f in flow:
        f.vars.remove("pos.x"); f.vars.remove("pos.y"); f.vars.remove("pos.z")
    nvar = len(flow[0].vars)
    for var in flow[0].vars:
        fp.write("%s\n" % var)
    fp.close()
    #
    # 3. Write .f file (function values)
    ffileName = rootName+(".t%04s" % tindx_str)+".f"
    ffileName = os.path.join(plotPath, ffileName)
    fp = open(ffileName, 'w')
    fp.write("%d\n" % len(flow))
    for f in flow:
        ni = f.ni; nj = f.nj; nk = f.nk
        if nk == 1:
            fp.write("%d %d %d\n" % (ni, nj, nvar))
        else:
            fp.write("%d %d %d %d\n" % (ni, nj, nk, nvar))
        for var in f.vars:
            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        fp.write("%20.12e\n" % uflowz(f.data[var][i,j,k]))
    #
    fp.close()
    return


#-----------------------------------------------------------------------------
# Profile extraction and writing.

def decode_range_from_string(range_str, min_value, max_value):
    """
    Returns the first and last indices corresponding to a range (with inclusive limits).
    """
    if range_str == ':':
        return min_value, max_value
    if range_str.find(':') >= 0:
        # We have a range specification to pull apart.
        items = range_str.split(':')
        first = int(items[0])
        if len(items) == 2:
            last = int(items[1])
        else:
            last = max_value
    else:
        # Presume that we have a single integer
        first = int(range_str)
        last = first
    #
    # negative indices start from the far end
    if first < 0: first = max_value + first + 1
    if first < min_value: first = min_value
    if first > max_value: first = max_value
    if last < 0: last = max_value + last + 1
    if last < min_value: last = min_value
    if last > max_value: last = max_value
    #
    return first, last

def write_profile_data(fileName, slice_list_str, tindx, nblock, grid, flow):
    """
    Write selected slices of cell data to a file in GnuPlot format.

    :param grid: list of StructuredGrid objects
    :param flow: list of StructuredGridFlow objects
    """
    fp = open(fileName, "w")
    flow[0].write_gnuplot_header(fp)
    slice_list = slice_list_str.split(';')
    # print "slice_list:", slice_list
    for slice_str in slice_list:
        bstr,istr,jstr,kstr = slice_str.split(',')
        bfirst,blast = decode_range_from_string(bstr, 0, nblock-1)
        for jb in range(bfirst,blast+1):
            kfirst,klast = decode_range_from_string(kstr, 0, flow[jb].nk-1)
            jfirst,jlast = decode_range_from_string(jstr, 0, flow[jb].nj-1)
            ifirst,ilast = decode_range_from_string(istr, 0, flow[jb].ni-1)
            print ("slice jb=%d i=%d:%d, j=%d:%d, k=%d:%d" %
                   (jb,ifirst,ilast,jfirst,jlast,kfirst,klast))
            for k in range(kfirst,klast+1):
                for j in range(jfirst,jlast+1):
                    for i in range(ifirst,ilast+1):
                        flow[jb].write_gnuplot_data_for_cell(fp, i, j, k)
    fp.close()
    return

def convert_string(slice_at_point_str, nblock, grid, flow):
    """
    Convert a slice-at-point string to a slice-list string
    that can be appended to another slice-at-list string.

    :param nblock: integer
    :param grid: list of StructuredGrid objects
    :param flow: list of StructuredGridFlow objects
    """
    slice_list_string = ""
    if len(slice_at_point_str) == 0: return slice_list_string
    slice_list = slice_at_point_str.split(';')
    for slice_str in slice_list:
        if len(slice_str) == 0: break
        bstr,index_pair,x,y,z = slice_str.split(',')
        bfirst,blast = decode_range_from_string(bstr, 0, nblock-1)
        pnt = Vector(float(x),float(y),float(z))
        for jb in range(bfirst,blast+1):
            ni = flow[jb].ni; nj = flow[jb].nj; nk = flow[jb].nk
            # Look for the nearest cell-centre in this block.
            # The search is along the index-direction not mentioned
            # in the index_pair string.
            if index_pair == 'ij' or index_pair == 'ji':
                i = 0; j = 0; kclose = None
                for k in range(nk):
                    for j in range(nj):
                        for i in range(ni):
                            x = flow[jb].data["pos.x"][i,j,k]
                            y = flow[jb].data["pos.y"][i,j,k]
                            z = flow[jb].data["pos.z"][i,j,k]
                            pcell = Vector(x,y,z)
                            d = vabs(pcell-pnt)
                            if kclose is None:
                                kclose = k; dclose = d
                            else:
                                if d < dclose:
                                    kclose = k; dclose = d;
                # We will have chosen one value for k,
                # whatever the position of pnt.
                slice_list_string += (";%d,0:%d,0:%d,%d" % (jb,ni,nj,kclose))
            elif index_pair == 'jk' or index_pair == 'kj':
                j = 0; k = 0; iclose = None
                for k in range(nk):
                    for j in range(nj):
                        for i in range(ni):
                            x = flow[jb].data["pos.x"][i,j,k]
                            y = flow[jb].data["pos.y"][i,j,k]
                            z = flow[jb].data["pos.z"][i,j,k]
                            pcell = Vector(x,y,z)
                            d = vabs(pcell-pnt)
                            if iclose is None:
                                iclose = i; dclose = d
                            else:
                                if d < dclose:
                                    iclose = i; dclose = d;
                # We will have chosen one value for i,
                # whatever the position of pnt.
                slice_list_string += (";%d,%d,0:%d,0:%d" % (jb,iclose,nj,nk))
            elif index_pair == 'ki'or index_pair == 'ik':
                i = 0; k = 0; jclose = None
                for k in range(nk):
                    for j in range(nj):
                        for i in range(ni):
                            x = flow[jb].data["pos.x"][i,j,k]
                            y = flow[jb].data["pos.y"][i,j,k]
                            z = flow[jb].data["pos.z"][i,j,k]
                            pcell = Vector(x,y,z)
                            d = vabs(pcell-pnt)
                            if jclose is None:
                                jclose = j; dclose = d
                            else:
                                if d < dclose:
                                    jclose = j; dclose = d;
                # We will have chosen one value for j,
                # whatever the position of pnt.
                slice_list_string += (";%d,0:%d,%d,0:%d" % (jb,ni,jclose,nk))
    return slice_list_string

def write_profile_along_line(fileName, slice_line_str, tindx, nblock, 
                             grid, flow, dimensions):
    """
    Write selected line of cell data to a file (in GnuPlot format).

    The line is defined as the straight line between p0 and p1, N sample points.
    There may be several lines specified, separated by semicolons. 
    """
    fp = open(fileName, "w")
    flow[0].write_gnuplot_header(fp)
    line_list = slice_line_str.split(';')
    for line_str in line_list:
        items = line_str.split(',')  # assume comma-separated values
        p0 = Vector(float(items[0]),float(items[1]),float(items[2]))
        p1 = Vector(float(items[3]),float(items[4]),float(items[5]))
        N = int(items[6])
        print "Profile along line p0=", p0, "p1=", p1, "N=", N
        ds = 1.0/float(N-1)
        jb=0; i=0; j=0; k=0  # initial guess for cell search
        for iN in range(N):
            s = iN*ds
            p = p0*(1.0-s) + p1*s
            print "    locate p=", p,
            jb, i, j, k = locate_cell_and_block(grid, flow, dimensions, i, j, k, jb, p.x, p.y, p.z)
            print "ijk=", i, j, k
            flow[jb].write_gnuplot_data_for_cell(fp, i, j, k)
    fp.close()
    return

#-----------------------------------------------------------------------------
# Comparing with a reference function.

def compute_difference_in_flow_data(f_ref, nblock, grid, flow, t):
    """
    Take difference with respect to a reference function.

    :param f_ref: function that returns flow data in a dictionary
    :param nblock: integer
    :param grid: list of StructuredGrid objects
    :param flow: list of StructuredGridFlow objects
    :param t:
    """
    for jb in range(nblock):
        ni = flow[jb].ni; nj = flow[jb].nj; nk = flow[jb].nk
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    x = flow[jb].data["pos.x"][i,j,k]
                    y = flow[jb].data["pos.y"][i,j,k]
                    z = flow[jb].data["pos.z"][i,j,k]
                    ref_flow_dict = f_ref(x, y, z, t)
                    for var in flow[jb].vars:
                        if var in ref_flow_dict.keys():
                            flow[jb].data[var][i,j,k] -= ref_flow_dict[var]
    return

def compute_difference_in_flow_data2(nblock, grid, flow, grid2, flow2, t):
    """
    Take difference with respect to a reference solution (grid2, flow2)

    Assumes same block topology (cell count, etc).
    """
    for jb in range(nblock):
        ni = flow[jb].ni; nj = flow[jb].nj; nk = flow[jb].nk
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    for var in flow[jb].vars:
                        if not(var in ["pos.x", "pos.y", "pos.z", "volume"]):
                            flow[jb].data[var][i,j,k] -= flow2[jb].data[var][i,j,k]
    return

def compute_volume_weighted_norms(nblock, grid, flow):
    """
    Returns a dictionary containing the integrated norms.

    This will make it easy to extract just the one or two of interest
    and might be useful when trying to automate convergence tests.
    """
    norms = {}
    norms["global"] = {}
    norms["per_block"] = {}
    for var in flow[0].vars:
        if not(var in ["pos.x", "pos.y", "pos.z", "volume"]):
            norms["global"]["%s"%var] = {}
            norms["global"]["%s"%var]["L1"] = 0.0
            norms["global"]["%s"%var]["L2"] = 0.0
            norms["global"]["%s"%var]["Linf"] = 0.0
            norms["global"]["%s"%var]["peak_pos"] = (0.0, 0.0, 0.0, 0)
    global_volume_sum = 0.0
    for jb in range(nblock):
        norms_for_block = {}
        ni = flow[jb].ni; nj = flow[jb].nj; nk = flow[jb].nk
        # First thing: compute block volume.
        block_volume_sum = 0.0
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    volume = flow[jb].data["volume"][i,j,k]
                    block_volume_sum += volume
                    global_volume_sum += volume
        # Second, compute norms for all flow variables.
        for var in flow[jb].vars:
            if not(var in ["pos.x", "pos.y", "pos.z", "volume"]):
                L1 = 0.0; L2 = 0.0; volume_sum = 0.0
                Linf= 0.0; peak_pos = (0.0, 0.0, 0.0)
                for k in range(nk):
                    for j in range(nj):
                        for i in range(ni):
                            x = flow[jb].data["pos.x"][i,j,k]
                            y = flow[jb].data["pos.y"][i,j,k]
                            z = flow[jb].data["pos.z"][i,j,k]
                            volume = flow[jb].data["volume"][i,j,k]
                            value = flow[jb].data[var][i,j,k]
                            L1 += volume * abs(value)
                            norms["global"]["%s"%var]["L1"] += volume*abs(value)
                            L2 += volume * value * value
                            norms["global"]["%s"%var]["L2"] += volume*value*value
                            if abs(value) > Linf:
                                Linf = abs(value)
                                peak_pos = (x, y, z)
                            if abs(value) > norms["global"]["%s"%var]["Linf"]:
                                norms["global"]["%s"%var]["Linf"] = abs(value)
                                norms["global"]["%s"%var]["peak_pos"] = (x, y, z, jb)
                norms_for_block["%s"%var] = {"L1":L1/block_volume_sum,
                                             "L2":math.sqrt(L2/block_volume_sum),
                                             "Linf":Linf,
                                             "peak_pos":peak_pos}
        norms_for_block["vol"] = block_volume_sum
        norms["per_block"]["block[%d]"%jb] = norms_for_block
    # Now, average to get the global norms for each variable.
    for var in flow[0].vars:
        if not(var in ["pos.x", "pos.y", "pos.z", "volume"]):
            norms["global"]["%s"%var]["L1"] = norms["global"]["%s"%var]["L1"]/global_volume_sum
            norms["global"]["%s"%var]["L2"] = math.sqrt(norms["global"]["%s"%var]["L2"]
                                                        /global_volume_sum)
    norms["global"]["vol"] = global_volume_sum
    return norms

def pretty_print_norms(norms, per_block_norm_list="", global_norm_list=""):
    """
    Print the norms to stdout so that we can easily find items in the output.

    Assume the particular structure of the dictionaries as defined in
    compute_volume_weighted_norms().

    The strings specifying which norms we want can specify several,
    each tuple separated by semicolons.
    To select a per-block norm, we specify a tuple "jb,var_name,norm_name".
    To select a global norm, we specify a tuple "var_name,norm_name".
    Available norms are "L1", "L2", "Linf" and "peak_pos".
    Of course, "peak_pos" goes with the "Linf" norm.

    These strings can also be "none" to indicate that no norms should be printed.
    """
    print "Volume-weighted norms, per-block:"
    if len(per_block_norm_list) == 0:
        # We've ask for norms without specifying which ones.
        # Print them all.
        for blk in norms["per_block"].keys():
            print "    ------- Begin", blk, "-------"
            for var in norms["per_block"][blk].keys():
                print "   ", var, norms["per_block"][blk][var]
            print "    ------- End block -------"
    elif per_block_norm_list.lower() == "none":
        print "    none requested"
    else:
        # Print selected norms.
        for item in per_block_norm_list.split(';'):
            jb, var_name, norm_name = item.split(',')
            block_name = "block["+jb+"]"
            print block_name, var_name, norm_name, \
                norms["per_block"][block_name][var_name][norm_name]
    print "Volume-weighted norms, global:"
    if len(global_norm_list) == 0:
        # We've ask for norms without specifying which ones.
        # Print them all.
        for var in norms["global"].keys():
            print "   ", var, norms["global"][var]
        print "------- End global norms ---------"
    elif global_norm_list.lower() == "none":
        print "    none requested"
    else:
        # Print selected global norms.
        for item in global_norm_list.split(';'):
            var_name, norm_name = item.split(',')
            print var_name, norm_name, norms["global"][var_name][norm_name]
    return
    
# -----------------------------------------------------------------------------
# radiation post-processing

def tangent_slab_along_slice(fileName, slice_list_str, tindx, nblock, grid, flow):
    """
    Perform a tangent-slab radiation calculation using the slices as the profile

    :param grid: list of StructuredGrid objects
    :param flow: list of StructuredGridFlow objects
    """
    
    cells = []; s = []
    slice_list = slice_list_str.split(';')
    # print "slice_list:", slice_list
    for slice_str in slice_list:
        bstr,istr,jstr,kstr = slice_str.split(',')
        bfirst,blast = decode_range_from_string(bstr, 0, nblock-1)
        for jb in range(bfirst,blast+1):
            kfirst,klast = decode_range_from_string(kstr, 0, flow[jb].nk-1)
            jfirst,jlast = decode_range_from_string(jstr, 0, flow[jb].nj-1)
            ifirst,ilast = decode_range_from_string(istr, 0, flow[jb].ni-1)
            print ("slice jb=%d i=%d:%d, j=%d:%d, k=%d:%d" %
                   (jb,ifirst,ilast,jfirst,jlast,kfirst,klast))
            for k in range(kfirst,klast+1):
                for j in range(jfirst,jlast+1):
                    for i in range(ifirst,ilast+1):
                    	cells.append( flow[jb].get_cell_data(i, j, k) )
                    	if len(cells)==1:
                    	    s.append( 0.0 )
                    	    p0 = Vector3( cells[-1]["pos.x"], cells[-1]["pos.y"], cells[-1]["pos.z"] )
                        else: 
                            p = Vector3( cells[-1]["pos.x"], cells[-1]["pos.y"], cells[-1]["pos.z"] )
                            s.append( vabs( p0 - p ) )
                            
    # tangent slab calculation
    nslabs = len(cells)
    rsm = create_radiation_spectral_model( "rad-model.lua" )
    gm = create_gas_model( "gas-model.lua" )
    nsp = gm.get_number_of_species()
    ntm = gm.get_number_of_modes()
    # T_i = cells[0]["T[0]"]
    # T_f = cells[-1]["T[0]"] + 0.5 * ( cells[-1]["T[0]"] - cells[-2]["T[0]"] )
    T_i = 0.0; T_f = 0.0
    # print "T_i = %e, T_f = %e\n" % ( T_i, T_f )
    TS = TS_data(rsm,nslabs,T_i,T_f)
    divq = []
    divq_OT = []
    print "Calculating the radiation spectra for all slabs..."
    b = Bar()
    b.len = 50
    rsm.prep_radiator_population_files()
    ds_vec = []
    # NOTE: must save gas-data structures to a list so that the pointers can 
    #       access the correct data
    Q_vec = []
    for islab in range(nslabs):
    	b.fill(int(float(islab)/float(nslabs)*100.0))
    	print b.show(), "\r",
    	sys.stdout.flush()
    	Q = Gas_data(gm)
    	Q.rho = cells[islab]["rho"]
    	for isp in range(nsp):
    	    sp = gm.species_name(isp)
    	    Q.massf[isp] = cells[islab]["massf[%d]-%s"%(isp,sp)]
    	for itm in range(ntm):
    	    Q.T[itm] = cells[islab]["T[%d]"%itm]
    	gm.eval_thermo_state_rhoT(Q)
    	divq.append( new_doublep() )
    	j_total = rsm.radiative_integrated_emission_for_gas_state(Q,True)
    	divq_OT.append( -j_total*4.0*math.pi )
    	if islab==0:
    	    ds = abs(s[1] - s[0])
    	else:
    	    ds = abs(s[islab] - s[islab-1])
    	TS.set_rad_point(islab,Q,divq[-1],s[islab],ds)
    	rsm.append_current_radiator_populations( s[islab] )
    	rsm.write_QSS_population_analysis_files( Q, islab )
    	ds_vec.append(ds)
    	Q_vec.append(Q)
    
    print "\nSolving the tangent-slab problem..."
    q_total = TS.solve_for_divq_OT()
    TS.F_.write_to_file("optically-thin-TS-flux-spectra.data")
    print "optically thin q_total = %0.3f W/cm2" % ( q_total * 1.0e-4 )
    q_total = TS.quick_solve_for_divq()
    TS.F_.write_to_file("approximate-TS-flux-spectra.data")
    print "approximate q_total = %0.3f W/cm2" % ( q_total * 1.0e-4 )
    q_total = TS.exact_solve_for_divq()
    TS.F_.write_to_file("exact-TS-flux-spectra.data")
    print "exact q_total = %0.3f W/cm2" % ( q_total * 1.0e-4 )
    ofile = open( fileName, "w" )
    ofile.write( "# Column 1: 1D cell location, s (m)\n" )
    ofile.write( "# Column 2: spatial step, ds (m)\n" )
    ofile.write( "# Column 3: Radiative divergence, divq (W/m**3)\n" )
    ofile.write( "# Column 4: Optically thin radiative divergence, divq (W/m**3)\n")
    for islab in range(nslabs):
        ofile.write("%e \t %e \t %e \t %e\n" % ( s[islab], ds_vec[islab], doublep_value(divq[islab]), divq_OT[islab] ) )
        delete_doublep(divq[islab])
    ofile.close()
    del rsm
    del gm

    return

