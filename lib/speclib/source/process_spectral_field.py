#!/usr/bin/env python
#    filename: process_spectral_field.py
# description: Summarize a spectral field spatially and spectrally

import sys
from spectral_field import *
from optparse import OptionParser
import pylab as plt
from code import interact

def usage():
    usage_str  =  "usage: process_spectral_field.py --spec-file=spec.txt\n"
    usage_str +=  "                                 --cal-file=cal.txt [none]\n"
    usage_str +=  "                                 --facility=X2/EAST [X2]\n"
    
    return usage_str
    
def main():
    # 0. Process command liine arguments    
    parser = OptionParser(usage=usage())
    parser.add_option( "--spec-file", action="store", type="string", dest="specFile",
                       help="Spectral field data file")
    parser.add_option( "--cal-file", action="store", type="string", dest="calFile",
                       help="Calibration data file (default to none)")
    parser.add_option( "--facility", action="store", type="string", dest="facility",
                       help="Facility the data is coming from (default to X2)")

    (options, args) = parser.parse_args()
         
    spec_data, IvL_all, IvX_all = process_spectral_field( options.specFile, options.calFile, options.facility )
    
    ans = raw_input("Would you like to enter interactive mode [y/n]? : ")
    if ans=="y":
        interact("Start interactive mode (Ctrl-D to return, or help(spec_data) for some hints)",local=locals())
    
    print "Done."

def process_spectral_field( specFile, calFile=None, facility=None ):
    # 0. Check inputs    
    if specFile == None:
        print usage(); sys.exit()

    # 1. Read in spectral field data
    if facility=="X2" or facility==None: spec_data = X2spectralField( specFile )
    elif facility=="EAST": spec_data = EASTspectralField( specFile )
    else:
        print "Facility: %s is not recognised"
        sys.exit()
    
    # 2. Apply calibration vector
    if calFile!=None: spec_data.apply_calibration_vector( calFile )
    else: spec_data.calibrated = True
    
    # 3. Present spatially averaged data
    IvL_all = spec_data.spatially_average()
    IvL_all.plot_data(xlabel="Wavelength, $\lambda$ (nm)",ylabel="Spectral radiance, $I_{\lambda}$ (W/cm2-um-sr)",label="IvL_all")
    
    # 3. Present spetrally integrated data
    IvX_all = spec_data.spectrally_integrate()
    IvX_all.plot_data(xlabel="Axial distance, x (m)",ylabel="Radiance, $I_{\lambda}$ (W/cm2-sr)",label="IvX_all")
    
    # 5. Return spec_data to the user
    return spec_data, IvL_all, IvX_all
    
if __name__=="__main__":
    main()
