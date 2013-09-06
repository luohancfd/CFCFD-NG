#!/bin/bash

#rad_model='QSS-nom-thin-DT.py'
rad_model='EQ-parade-DT.py'

# 1. Create the radiation model
radmodel.py -i $rad_model -L rad-model.lua

# 2. Run e3prep
e3prep.py --job=hayabusa

