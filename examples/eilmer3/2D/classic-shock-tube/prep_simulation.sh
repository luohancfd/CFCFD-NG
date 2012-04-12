#/bin/bash

if [ -f ./cea-lut-air.lua.gz ]
then
    echo "Found LUT file already in place."
else
    echo "Generate LUT file for air."
    build-cea-lut.py --gas=air
fi
e3prep.py --job=cst --do-svg
