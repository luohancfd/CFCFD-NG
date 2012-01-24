#! /bin/sh
# file: setup_LUT_plus_composite_ideal_gas.sh
#
echo "Build look-up table for air."
build-cea-lut --case=air
echo "We should now have a look-up table for air"

# Note that, with changes to the functions in l_script.py, the following  actions
# are now done by the Python code when the user script slects the gas model.

echo "Set up composite ideal gas."
echo "model = 'ideal gas'\nspecies = {'Ar','He','N2','air'}" > gas-tmp.inp
gasfile gas-tmp.inp gas-model-tmp-0.lua

echo "Add look-up table file to gas-model file."
sed "1alut_file = 'cea-lut-air.lua.gz'" gas-model-tmp-0.lua > gas-model-tmp-1.lua
sed "s/composite gas/LUT-plus-composite/" gas-model-tmp-1.lua > gas-model.lua
echo "We should now have gas-model.lua file ready for simulation."

