#!/bin/bash

python make_source_terms.py
e4shared --prep --job=mms
e4shared --run --job=mms --max-cpus=1
e3post.py --job=mms --tindx=20 --ref-function=analytic_solution.py  --global-norm-list="rho,L2"

