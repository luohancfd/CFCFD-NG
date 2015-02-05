#!/bin/bash
# test.sh
e3prep.py --job=test --json
cp ./sample-data/ideal-air-gas-model.lua ./gas-model.lua
./e4shared --run --job=test --verbosity=2
