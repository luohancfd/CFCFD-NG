#!/usr/bin/env python
## \file radmodel.py
## \ingroup radiation
## \brief Python program to create a radiation input file
##
## \author Daniel F Potter
## \version 13-Jan-2013: reborn as radmodel.py (previously script_rad2.py)

import os
import sys
from datetime import datetime
from getopt import getopt

import radiator_library as rl
import parade_radiator_library as prl

tab = rl.tab

sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory

class GlobalRadData(object):
    """Python class to organize the global data for a radiation system.

    The user's script does not create this object but rather just alters
    attributes of the global object.
    """
    count = 0
    def __init__(self):
        self.spectral_model = "none"
        self.radiators = []
        self.lambda_min = 0.0
        self.lambda_max = 0.0
        self.spectral_points = 0
        self.spectral_blocks = 1
        self.adaptive_spectral_grid = False
        self.transport_model = "none"
        self.nrays = 0
        self.spectrally_resolved = True
        self.upper_escape_factor = 0.0
        self.lower_escape_factor = 0.0
        self.optical_switch = 0.0
        self.electronic_mode_factor = 0.0
        self.absorption = "partitioned energy"
        self.binning = "none"
        self.N_bins = 0
        self.exact_formulation = False

        return
    
    def write_LUA_file( self, ofile_str, ifile_str ):
        # Firstly check that there are radiators present
        if ( len(self.radiators)==0 ):
            print "No radiators have been requested!"
            sys.exit()
        ofile = open( ofile_str, "w" )
        # file header
        now = datetime.now()
        ofile.write("-- File: %s\n" % ( ofile_str ) )
        ofile.write("-- This file was automatically created by 'script_rad2.py' from \n")
        ofile.write("-- input script '%s' at %s\n" % ( ifile_str, str(now) ) )
        ofile.write("\n")
        # spectral data
        ofile.write("spectral_data = {\n")
        ofile.write(tab+"spectral_model = '%s',\n" % self.spectral_model )
        if self.spectral_model=="parade":
            ofile.write(tab+"control_template = 'control.template',\n" )
            snbopt_template = 'none'
            for rad in self.radiators:
                if rad.type=="triatomic_radiator" and rad.band_model=="SNB":
                    snbopt_template = "SNBOPT.template"
            ofile.write(tab+"snbopt_template = '%s',\n" % snbopt_template )
            iT = self.radiators[0].iT; iTe = self.radiators[0].iTe
            for rad in self.radiators[1:]:
                if rad.iT!=iT:
                    print "iT not consistent amongst radiators"
                    sys.exit()
                if rad.iTe!=iTe:
                    print "iTe not consistent amongst radiators"
                    sys.exit()
            
            ofile.write(tab+"iT = %d,\n" % iT )
            ofile.write(tab+"iTe = %d,\n" % iTe )
        ofile.write(tab+"radiators = { ")
        for rad in self.radiators:
            ofile.write("'%s', " % rad.name )
        ofile.write("},\n")
        if self.lambda_max < self.lambda_min:
            print "lambda_min is greater than lambda_max!"
            sys.exit()
        ofile.write(tab+"lambda_min = %f,\n" % self.lambda_min )
        ofile.write(tab+"lambda_max = %f,\n" % self.lambda_max ) 
        if self.spectral_points<1:
            print "spectral_points is less than 1!"
            sys.exit()
        ofile.write(tab+"spectral_points = %d,\n" % self.spectral_points )
        if self.spectral_blocks<1:
            print "spectral_blocks is less than 1!"
            sys.exit()
        ofile.write(tab+"spectral_blocks = %d,\n" % self.spectral_blocks )
        ofile.write(tab+"adaptive_spectral_grid = %s,\n" % str(self.adaptive_spectral_grid).lower() )
        ofile.write("}\n\n")
        # transport data
        ofile.write("transport_data = {\n")
        ofile.write(tab+"transport_model = '%s',\n" % self.transport_model )
        ofile.write(tab+"spectrally_resolved = %s,\n" % str(self.spectrally_resolved).lower() )
        ofile.write(tab+"electronic_mode_factor = %d,\n" % self.electronic_mode_factor )
        if self.transport_model=="discrete transfer" or self.transport_model=="monte carlo":
            if self.nrays<1:
                print "nrays is less than 1!"
                sys.exit()
            ofile.write(tab+"nrays = %d,\n" % self.nrays )
            ofile.write(tab+"clustering = 'none',\n" )
        if self.transport_model=="discrete transfer":
            ofile.write(tab+"binning = '%s',\n" % self.binning )
            ofile.write(tab+"N_bins = '%d',\n" % self.N_bins )
        elif self.transport_model=="monte carlo":
            ofile.write(tab+"absorption = '%s',\n" % self.absorption )
        elif self.transport_model=="optically variable":
            ofile.write(tab+"optical_switch = %f,\n" % self.optical_switch )
            ofile.write(tab+"lower_escape_factor = %f,\n" % self.lower_escape_factor )
            ofile.write(tab+"upper_escape_factor = %f,\n" % self.upper_escape_factor )
        elif self.transport_model=="tangent slab":
            ofile.write(tab+"exact_formulation = %s,\n" % str(self.exact_formulation).lower() )
        ofile.write("}\n\n")
            # radiator data
        for rad in self.radiators:
            ofile.write("%s\n" % rad.get_LUA_string() )
        # finished
        ofile.close()
    
    def request_radiator( self, rrad ):
        if self.spectral_model=="photaura":
            available_radiators = rl.available_radiators
        elif self.spectral_model=="parade":
            available_radiators = prl.available_radiators
        else:
            print "Spectral model: %s not understood" % self.spectral_model
            sys.exit()
        for arad in available_radiators.keys():
            if rrad==arad: 
                self.radiators.append( available_radiators[arad] )
                return self.radiators[-1]
        
        print "Requested radiator: %s was not found.\n" % rrad
        print "Available radiators are: ", available_radiators.keys()
        sys.exit()

    def create_parade_template_files( self, data_path ):
        ffrad = False; fbrad = False; alrad = False; mbrad = False; snb = False
        for rad in self.radiators:
            if rad.type=="electron_radiator":
                ffrad = True; fbrad = True
            elif rad.type=="atomic_radiator":
                alrad = True
            elif rad.type=="diatomic_radiator" or rad.type=="triatomic_radiator":
                mbrad = True
            if rad.type=="triatomic_radiator" and rad.band_model=="SNB":
                snb = True
            ofile = open("control.template","w")
            ofile.write("c       This PARADE 3.1 control template file was automatically created\n")
            ofile.write("c       by TRT_tools.py at %s\n" % datetime.now() )    
            ofile.write("c       1. Spectrum control data:\n")
            ofile.write("c\n")
            ofile.write(" %d\t\twavlo   [A]\n" % (self.lambda_min*10) )
            ofile.write(" %d\t\twavhi   [A]\n" % (self.lambda_max*10) )
            ofile.write(" %d\t\tnpoints [A]\n" % self.spectral_points   )
            ofile.write("c\n")
            ofile.write("c       2. Parameters for adaptive wavelength discretisation (not used):\n")
            ofile.write("  0.0005       minimum distance between adjacent points\n")
            ofile.write("   100.        integration limit for line shape\n")
            ofile.write("  0.001        minimum fraction for use of energy level\n")
            ofile.write("c       3. Switches for radiation mechanisms:\n")
            ofile.write("  9            number of control switches\n")
            if ffrad:
                ofile.write("  'Y'          free-free radiation (y/n)\n")
            else:
                ofile.write("  'N'          free-free radiation (y/n)\n")
            if fbrad:
                ofile.write("  'Y'          free-bound radiation (y/n)\n")
            else:
                ofile.write("  'N'          free-bound radiation (y/n)\n")
            if alrad:
                ofile.write("  'Y'          atomic line radiation (y/n)\n")
            else:
                ofile.write("  'N'          atomic line radiation (y/n)\n")
            if mbrad:
                ofile.write("  'Y'          molecular band radiation (y/n)\n")
            else:
                ofile.write("  'N'          molecular band radiation (y/n)\n")
            ofile.write("  'N'          read 'parade.rad' if available (y/n)\n")
            ofile.write("  'N'          adaptive wavelength discretisation (y/n)\n")
            ofile.write("  'N'          equal wavelength increments (y/n)\n")
            ofile.write("  'Y'          equal frequency increments (y/n)\n")
            if snb:
                ofile.write("  'Y'          Statistical Narrow-Band randiation (y/n)\n")
            else:
                ofile.write("  'N'          Statistical Narrow-Band randiation (y/n)\n")
            ofile.write("c 3(bis). Switches for output options\n")
            ofile.write("  1         iout (row number for main outputs)\n")
            ofile.write("  1        jout (column number for main outputs)\n")
            ofile.write("  3            number of output switches\n")
            ofile.write("  'N'          'par_res.imo' for each cell (y/n)\n")
            ofile.write("  'N'          time integrated emission coefficient until cell number\n")
            ofile.write("  2            debug level (0: minimum output, 2: maximum output)\n")
            ofile.write("  60           max cell number for time integration (shock tube test)\n")
            ofile.write("c\n")
            ofile.write("c       4. rad(y/n) remark  at. spec  tt tr  tv  te    rad.file\n")
            mol_index = 0
            for rad in self.radiators:
                if rad.atoms>1: mol_index+=1
                ofile.write(rad.get_parade_string(mol_index,data_path))
            ofile.close()
            
        # may also need an SNBOPT template file
        if snb:
            ofile = open("SNBOPT.template","w")
            ofile.write("4\n")
            ofile.write("9       17      29      70\n")
            ofile.write("x       i\n")
            ofile.close()

def main():
    from optparse import OptionParser

    usage =  "usage: radmodel.py -i rad_desc.py|--input-script=rad_desc.py\n"
    usage += "                   -L LUA_output.lua|--LUA-file=LUA_output.lua"
    parser = OptionParser(usage=usage)
    parser.add_option( "-i", "--input-script",
                       action="store", type="string", dest="inFile",
                       help="input Python script for radiation description")
    parser.add_option( "-L", "--LUA-file",
                       action="store", type="string", dest="LUAFile",
                       help="output configuration file for 'librad' C++ module in LUA format")

    (options, args) = parser.parse_args()
    
    if options.inFile==None:
        print usage
        sys.exit()
    
    gdata = GlobalRadData()
    execfile(options.inFile)
    
    if ( options.LUAFile != None ):
        gdata.write_LUA_file( options.LUAFile, options.inFile )

    
if __name__ == '__main__':
    main()

