#!/bin/bash

# 1. Create the radiation model
script_rad2.py -i QSS-nom-thin-DT.py -L rad-model.lua

# 2. Run e3prep
e3prep.py --job=hayabusa

