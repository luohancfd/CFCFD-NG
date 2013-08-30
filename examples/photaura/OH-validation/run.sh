#!/bin/bash

# 1. Setup the radiation model LUA input file
radmodel.py -i OH-radiators.py -L rad-model.lua

# 2. Run EQ_spectra.py
EQ_spectra.py --input-file=OH.py

# The files coefficient_spectra.txt, coefficient_spectra_with_AF.txt
# and intensity_spectra.txt should now be in the working directory.

# 3. Plot the results

# 3a. photaura vs lifbase, non-dimensional comparison (normalised to peaks) 
gnuplot plot-Photaura-LB.gnu
epstopdf emission-spectra-Photaura-LB.eps

# 3b. photaura vs spartan, dimensional comparison 
gnuplot plot-Photaura-SPARTAN.gnu
epstopdf emission-Photaura-SPARTAN.eps
epstopdf absorption-Photaura-SPARTAN.eps

# 3c. photaura vs lifbase vs spartan, normalised to peaks 
gnuplot plot-spectra.gnu
epstopdf emission-spectra.eps

# 3d. photaura vs experiment
gnuplot plot-Photaura-exp.gnu
epstopdf intensity-Photaura-exp.eps
epstopdf intensity-Photaura-exp-SL1.eps
epstopdf intensity-Photaura-exp-SL2.eps
epstopdf intensity-Photaura-exp-SL3.eps

evince *pdf
