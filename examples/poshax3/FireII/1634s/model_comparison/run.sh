#!/bin/bash

# 1. Two-temperature gas-model calculations
cd 2T/

# 1.1 Park chemical kinetics
cd Park/
./run.sh
cd ../

cd ../

# 2. Three-temperature gas-model calculations
cd 3T/

# 2.1 Park chemical kinetics
cd Park/
./run.sh
cd ../

# 2.1 Park chemical kinetics with Knab CVCV model
cd Knab/
./run.sh
cd ../

cd ../

# 3. Plot comparison of electron number densities
gnuplot profiles.gplot
