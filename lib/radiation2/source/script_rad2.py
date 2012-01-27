#!/usr/bin/env python
## \file script_rad2.py
## \ingroup radiation2
## \brief Python program to create a radiation2 input file
##
## \author Daniel F Potter
## \version 13-Aug-2009

import os
import sys
from datetime import datetime
from getopt import getopt

import radiator_library2 as rl

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

        if GlobalRadData.count >= 1:
            raise Exception, "Already have a GlobalNoneqData object defined."
	self.spectral_model = "none"
	self.radiators = []
	self.lambda_min = 0.0
	self.lambda_max = 0.0
	self.spectral_points = 0
	self.spectral_blocks = 1
	self.transport_model = "none"
	self.nrays = 0
	self.spectrally_resolved = True
	self.upper_escape_factor = 0.0
	self.lower_escape_factor = 0.0
	self.optical_switch = 0.0

        GlobalRadData.count = 1
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
	ofile.write("}\n\n")
	# transport data
	ofile.write("transport_data = {\n")
	ofile.write(tab+"transport_model = '%s',\n" % self.transport_model )
	ofile.write(tab+"spectrally_resolved = %d,\n" % self.spectrally_resolved )
	if self.transport_model=="discrete transfer" or self.transport_model=="monte carlo":
	    if self.nrays<1:
	        print "nrays is less than 1!"
		sys.exit()
	    ofile.write(tab+"nrays = %d,\n" % self.nrays )
	    ofile.write(tab+"clustering = 'none',\n" )
            ofile.write(tab+"binning = 'none',\n" )
	elif self.transport_model=="optically variable":
	    ofile.write(tab+"optical_switch = %f,\n" % self.optical_switch )
	    ofile.write(tab+"lower_escape_factor = %f,\n" % self.lower_escape_factor )
	    ofile.write(tab+"upper_escape_factor = %f,\n" % self.upper_escape_factor )
	ofile.write("}\n\n")
        # radiator data
	for rad in self.radiators:
	    ofile.write("%s\n" % rad.get_LUA_string() )
	# finished
	ofile.close()
	
    def request_radiator( self, rrad ):
	for arad in rl.available_radiators.keys():
	    if rrad==arad: 
		self.radiators.append( rl.available_radiators[arad] )
		return self.radiators[-1]
		
        print "Requested radiator: %s was not found.\n" % rrad
        print "Available radiators are: ", rl.available_radiators.keys()
        sys.exit()
		
		
gdata = GlobalRadData()

def main():
    from optparse import OptionParser

    usage =  "usage: script_rad2.py -i rad_desc.py|--input-script=rad_desc.py\n"
    usage += "                      -L LUA_output.lua|--LUA-file=LUA_output.lua"
    parser = OptionParser(usage=usage)
    parser.add_option( "-i", "--input-script",
                       action="store", type="string", dest="inFile",
                       help="input Python script for radiation description")
    parser.add_option( "-L", "--LUA-file",
                       action="store", type="string", dest="LUAFile",
                       help="output configuration file for 'librad2' C++ module in LUA format")

    (options, args) = parser.parse_args()
    
    if options.inFile==None:
        print usage
	sys.exit()
    
    execfile(options.inFile)
	
    if ( options.LUAFile != None ):
        gdata.write_LUA_file( options.LUAFile, options.inFile )

    
if __name__ == '__main__':
    main()

