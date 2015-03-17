#!/bin/bash

python make_source_terms.py
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=mms
e4shared --run --job=mms
e3post.py --job=mms --tindx=20 --ref-function=analytic_solution.py  --global-norm-list="rho,L2"

