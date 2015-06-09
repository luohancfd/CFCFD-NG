#!/bin/bash
python make_source_terms.py
e4shared --job=mms --prep
e4shared --job=mms --run
e3post.py --job=mms --tindx=20 --ref-function=analytic_solution.py  --global-norm-list="rho,L2;rho,Linf" | tail -2 > rho-norms.txt
