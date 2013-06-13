#!/usr/bin/env python
import EQ_spectra
import radmodel
import os
from getopt import getopt, GetoptError
import sys

longOptions = ["help", "with-photaura", "with-spradian", "with-parade", "with-equilibrium-air" ]

def printUsage():
    print ""
    print "Usage: test-spectral-models.py [--help] [--with-photaura] [--with-equilibrium-air] [--with-spradian] [--with-parade]"
    print "e.g. test-spectral-models.py --with-photaura --with-equilibrium-air --with-spradian"
    print ""
    return

try:
    userOptions = getopt(sys.argv[1:], [], longOptions)
except GetoptError, e:
    print "One (or more) of your command-line options was no good."
    print "    ", e
    printUsage()
    sys.exit(1)
uoDict = dict(userOptions[0])
if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
    printUsage()
    sys.exit(0)
#
with_photaura = uoDict.has_key("--with-photaura")
#
with_equilibrium_air = uoDict.has_key("--with-equilibrium-air")
#
with_spradian = uoDict.has_key("--with-spradian")
#
with_parade = uoDict.has_key("--with-parade")
           
models = dict()
if with_photaura:
    models["photaura"] = "photaura-air-EQ-radiators.py"
if with_equilibrium_air:
    models["equilibrium air"] = "equilibrium-air.lua"
if with_spradian:
    models["spradian"] = "spradian-air-EQ-radiators.lua"
if with_parade:
    models["parade"] = "parade-air-EQ-radiators.py"
    
if len(models.keys())==0:
    print "No models requested for testing!"
    printUsage()
    sys.exit()
           
intensity_results = dict()
expected_intensities = { "photaura"        : 1.249113903993e+05,
                         "equilibrium air" : 1.076e+05,
                         "spradian"        : 2.347982840522e+05,
                         "parade"          : 8.481522557937e+04  }
time_results = dict()
           
for model in models.keys():
    # 1.Create the lua file for this model
    gdata = radmodel.GlobalRadData()
    ifile = models[model]
    if ifile.find(".py")>=0:
        execfile(ifile)
        gdata.write_LUA_file( "rad-model.lua", ifile )
    else:
        os.system("cp %s rad-model.lua" % ifile)
        
    # 2. Run EQ-spectra
    input_data = EQ_spectra.EqSpectraInputData()
    input_data.rad_model_file = "rad-model.lua"
    input_data.species_list = [ "N2", "N2_plus", "NO", "NO_plus", "O2", "O2_plus", "N", "N_plus", "O", "O_plus", "e_minus" ]
    input_data.mass_fractions = [ 0.767, 0.0, 0.0, 0.0, 0.233, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    input_data.shock_speed = 10000.0 # m/s
    input_data.gas_pressure = 0.1*133.33 # Pa
    input_data.gas_temperature = 296 # K
    input_data.path_length = 10.0e-2 # m
    input_data.apparatus_fn = "none"
    input_data.planck_spectrum = False
    input_data.show_plots = False
    I_total, delta_t = EQ_spectra.run_calculation(input_data)
    intensity_results[model] = I_total
    time_results[model] = delta_t
    
def percent_diff(x, y):
    return 100.0 * abs(x-y) / max(abs(x),abs(y))
    
print "\n-------------------------------------------------------------------------------------------"
ntests = 0
nfail = 0
npass = 0
for model in models.keys():
    ntests += 1
    error = percent_diff(intensity_results[model],expected_intensities[model])
    if error >= 0.1:
        flag = "FAILED"
        nfail += 1
    else:
        flag = "PASSED"
        npass += 1
    print "model: %s \t I_total: %0.2f W/cm2-sr \t wall-time: %0.2f s \t [%s]" % ( model, intensity_results[model]*1.0e-4, time_results[model], flag )
print "-------------------------------------------------------------------------------------------"
print "Total tests: %d \t passed: %d \t failed: %d" % ( ntests, npass, nfail )
print "-------------------------------------------\n"
print "Test complete."

