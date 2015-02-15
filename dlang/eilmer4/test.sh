#!/bin/bash
# test.sh
# We need to shuffle the gas-model files because Eilmer4 uses
# Rowan's new gas model format while Eilmer3 uses the older format.

e3prep.py --job=test --json
mv ./gas-model.lua ./gas-model-save.lua

cp ./sample-data/ideal-air-gas-model.lua ./gas-model.lua
./e4shared --run --job=test --verbosity=1

cp ./gas-model-save.lua ./gas-model.lua
e3post.py --job=test --tindx=all --vtk-xml
