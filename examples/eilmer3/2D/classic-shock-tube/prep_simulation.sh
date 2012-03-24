#/bin/bash

if [ -f ./cea-lut-air.lua.gz ]
then
    echo "Found LUT file already in place."
else
    echo "Generate LUT file for air."
    build-cea-lut --case=air
fi
e3prep.py --job=cst
