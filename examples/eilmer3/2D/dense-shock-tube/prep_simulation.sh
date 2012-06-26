#/bin/bash

if [ -f ./REFPROP-lut-CO2.FLD.lua.gz ]
then
    echo "Found LUT file already in place."
else
    echo "Generate LUT file for air."
    build-REFPROP-lut.py --fluid=CO2.FLD --bounds="300,1000,-3,1.54"
fi
e3prep.py --job=dst --do-svg
