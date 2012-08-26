#! /bin/sh
# here we test the Noneq_thermal_behaviour class by:
# - setting up a 2T 5sp air gas-model 
# - evaluating thermo properties over a wide temperature
#   range setting Tv=T (thermal equilibrium)
# - comparing results with that of cea2 and a 
#   Perfect_thermal_behaviour model
./clean.sh
noneq-test.py
gnuplot plot.gnu
