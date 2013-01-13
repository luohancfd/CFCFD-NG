#!/bin/bash

# 1. Create the radiation model
radmodel.py -i QSS-nom-thin-DT.py -L rad-model.lua

# 2. Run e3prep
e3prep.py --job=hayabusa

