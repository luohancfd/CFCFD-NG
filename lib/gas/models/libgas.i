/** \file libgas.i
 *  \ingroup libgas 
 *  \brief SWIG interface header file for C++ gas model library.
 *  \author PJ, RJG
 *  \version 07-Feb-2006
 *  \version 14-May-2008 moved to its new home lib/gas/models/
 */
%define DOCSTRING
"Python interface to the gas model classes, now implemented in C++"
%enddef
%module(docstring=DOCSTRING) libgas

%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"

// The following magic allows us to pass a list of numbers
// to be collected as a C++ vector in the CFlowCondition constructor.
#ifndef STD_VECTOR_TEMPLATES_ALREADY_DEFINED
%template(vectori) std::vector<int>;
%template(vectord) std::vector<double>;
#define STD_VECTOR_TEMPLATES_ALREADY_DEFINED
#endif
%template(map_str_int) std::map<std::string, int>;
%{
#include <string>
#include <vector>
#include <map>
#include <cstdio>
#include "gas_data.hh"
#include "gas-model.hh"
#include "../../util/source/lua_service.hh"
#include "../../nm/source/Richardson_extrapolation.hh"
#include "user-defined-gas-model.hh"
#include "lservice_gas_data.hh"
#include "physical_constants.hh"
#include "../kinetics/reaction-rate-coeff.hh"
#include "../kinetics/reaction-update.hh"
#include "../kinetics/generalised-Arrhenius.hh"
#include "../kinetics/Park-nonequilibrium.hh"
#include "../kinetics/MarroneTreanor-dissociation.hh"
#include "../kinetics/Macheret-dissociation.hh"
#include "../kinetics/Knab-molecular-reaction.hh"
%}

%ignore Gas_data::operator=;
%include "gas_data.hh"
%include "gas-model.hh"
%include "physical_constants.hh"
%include "../kinetics/reaction-rate-coeff.hh"
%include "../kinetics/reaction-update.hh"
%include "../kinetics/generalised-Arrhenius.hh"
%include "../kinetics/Park-nonequilibrium.hh"
%include "../kinetics/MarroneTreanor-dissociation.hh"
%include "../kinetics/Macheret-dissociation.hh"
%include "../kinetics/Knab-molecular-reaction.hh"


%pythoncode%{
import sys
import os
import re
from numpy import zeros

R_u = PC_R_u
R_u_kmol = PC_R_u_kmol
Avogadro = PC_Avogadro
kB = PC_k_SI
P_atm = PC_P_atm
p_atm = PC_P_atm
T_atm = PC_T_ref

def set_molef(Q, g, X):
    """
    Set mole-fractions vector from a dictionary or list.
    """ 
    nsp = g.get_number_of_species()
    if type(X) is dict:
        molef = zeros(nsp)
        for species in X.keys():
            species_exists = False
            for isp in range(nsp):
                if species == g.species_name(isp):
                    molef[isp] = X[species]
                    species_exists = True
            if (not species_exists):
                print "input error: species %s does not exist." % species
                sys.exit(1)
    elif type(X) is list:
        molef = []
        for isp in range(nsp):
            molef.append(X[isp])
    else:
        print "input error: species should be provided in a dict or list."
        sys.exit(1)
    # convert to massf
    massf = convert_molef2massf(molef, g.M())
    for isp in range(nsp):
	Q.massf[isp] = massf[isp]
    return

def set_massf(Q, g, mf):
    """
    Set mass fraction vector from a dictionary or list.
    """
    nsp = g.get_number_of_species()
    if type(mf) is dict:
        massf = zeros(nsp)
        for species in mf.keys():
            species_exists = False
            for isp in range(nsp):
                if species == g.species_name(isp):
                    massf[isp] = mf[species]
                    species_exists = True
            if (not species_exists):
                print "input error: species %s does not exist." % species
                sys.exit(1)
    elif type(mf) is list:
        massf = []
        for isp in range(nsp):
            massf.append(mf[isp])
    else:
        print "input error: species should be provided in a dict or list."
        sys.exit(1)
    # normalise
    sum_mf = sum(massf)
    for isp in range(nsp):
	Q.massf[isp] = massf[isp]/sum_mf
    return

# Functions moved from e3_gas.py
# Author: Rowan J. Gollan
# Date: 28-July-2008
# Moved here, simplified and updated for LUT by PJ, May 2010.


def create_gas_file(model, species, fname="gas-model.lua", lut_file=None):
    """
    Write the input file for Rowan's gasfile.lua program.

    That program is then invoked to look up the species database
    and put together the detailed gas-model.lua file that is 
    input to the main C++ codes.

    Return a pointer to a gas model object which may be
    optionally caught by the caller.

    Input:
    model   : (string) name of the gas model as specified in gasfile.lua.
    species : list of species names (strings)
    fname   : (string) name of the file to be output be gasfile.lua
    lut_file: (string) name of the Look-up table file, if any.

    Note that the species list does not contain the LUT species
    which will the automatically prepended if the lut_fule is specified.
    This addition only works if the original gas model is a composite-gas.
    The resulting gas model is LUT-plus-composite.
    """
    # Prepare the basic gas model file by invoking Rowan's gasfile.lua program.
    tmpfile = "gas-tmp.txt"
    fp = open(tmpfile, "w")
    fp.write("model = '%s'\n" % model)
    fp.write("species = {")
    for spec in species:
        fp.write("'%s', " % spec)
    fp.write("}\n")
    fp.close()
    cmd = "%s %s %s" % ("gasfile", tmpfile, fname)
    os.system(cmd)
    cmd = "rm " + tmpfile
    os.system(cmd)
    if lut_file:
        # We need to modify the generated gas-model.lua file to include LUT details.
        fp = open(fname, "r")
        lines_of_text = fp.readlines()
        fp.close()
        if " ".join(lines_of_text).find("composite gas") < 0:
            print "A LUT gas can only be added to a composite-gas model."
            sys.exit(-1) 
        lines_of_text.insert(1, "-- modified by create_gas_file() to add LUT\n")
        lines_of_text.insert(2, "lut_file = '%s'\n" % lut_file)
        fp = open(fname, "w")
        for line in lines_of_text:
            if line.find("composite gas") >= 0:
                line = line.replace("composite gas", "LUT-plus-composite")
            fp.write(line)
        fp.close()
    #
    return create_gas_model(fname)


def change_ideal_gas_attribute(species_name, attribute_name, new_value,
                               fname="gas-model.lua"):
    """
    A convenient way to change gas attributes from the user-input script.

    This function assumes that the user has selected a simple
    ideal gas, and as such, the gas configuration file is gas-model.lua.
    So first it tests for the existence of this file, and then performs
    the appropriate substitution.
    """
    if not os.path.isfile(fname):
        print "change_ideal_gas_attribute()"
        print "This function assumes that a file", fname, "is present."
        print fname, "has NOT been found."
        print "Bailing out!"
        sys.exit(1)

    import shutil
    shutil.move(fname, fname+"~")

    destination = open(fname, "w")
    source = open(fname+"~", "r")

    pat = species_name + '.' + attribute_name
    pat = re.compile(pat)
    val_str = re.compile(' *value')

    # Loop over lines in file looking for
    # attribute to amend
    found = False
    while 1:
        line = source.readline()
        if not line:
            break
        destination.write(line)
        if pat.match(line):
            found = True
            while 1:
                line = source.readline()
                if not line:
                    print "There is a problem finding the value string for attribute:"
                    print species_name + '.' + attribute_name
                    print "Leaving", fname, "unchanged."
                    source.close()
                    destination.close()
                    shutil.move(fname+"~", fname)
                    print "Bailing out!"
                    sys.exit(1)

                if val_str.match(line):
                    # Do the substituion.
                    destination.write("  value = %e,\n" % new_value)
                    break
                else:
                    destination.write(line)
    source.close()
    destination.close()
    if not found:
        print "The requested attribute was not found in", fname, ":"
        print species_name + '.' + attribute_name
        print "Leaving", fname, "unchanged."
        shutil.move(fname+"~", fname)
    os.remove(fname+"~")
    return


# Some methods added to the Gas_model class for the convenience of the user.
# PJ, 01-Feb-2011

def get_species_names(self):
    nsp = self.get_number_of_species()
    return [self.species_name(isp) for isp in range(nsp)]
Gas_model.get_species_names = get_species_names

def to_list(self, value_dict):
    """
    Returns a filled-in list of values for all of the species.
    
    The values may be mass-fractions or mole-fractions.
    self: Gas_model
    value_dict: the dictionary containing the values, with species names as keys.
    """
    nsp = self.get_number_of_species()
    value_list = [0.0,] * nsp
    species_names = self.get_species_names()
    for species in value_dict.keys():
        value_list[species_names.index(species)] = float(value_dict[species])
    return value_list
Gas_model.to_list = to_list

def to_molef(self, mf):
    """
    Converts from mass-fractions to mole-fractions.

    self : Gas_model
    mf   : mass fractions in dictionary or list format
    """
    if type(mf) is dict:
        mf = self.to_list(mf)
    return convert_massf2molef(mf, self.M())
Gas_model.to_molef = to_molef

def to_massf(self, mf):
    """
    Converts from mole-fractions to mass-fractions.

    self : Gas_model
    mf   : mole fractions in dictionary or list format
    """
    if type(mf) is dict:
        mf = self.to_list(mf)
    return convert_molef2massf(mf, self.M())
Gas_model.to_massf = to_massf

def set_massf_vector(self, Q, mf):
    """
    Set the full mass fractions vector.
    """
    set_massf(Q, self, mf)
    return
Gas_model.set_massf = set_massf_vector

def set_molef_vector(self, Q, mf):
    """
    Set the full mass fractions vector from mole fractions.
    """
    set_molef(Q, self, mf)
    return
Gas_model.set_molef = set_molef_vector

def get_atomic_constituents(self, isp):
    """
    Return the atomic constituents as a dict.
    """
    map = map_str_int({})
    self.atomic_constituents(isp, map)
    return dict(map)
Gas_model.get_atomic_constituents = get_atomic_constituents

def list_all_atoms(self):
    """
    Return a list of all atomic species in mix.
    """
    nsp = self.get_number_of_species()
    atoms = set({})
    for isp in range(nsp):
        d = self.get_atomic_constituents(isp)
	for k in d.keys():
            atoms.add(k)
    return list(atoms)
Gas_model.list_all_atoms = list_all_atoms


%}
