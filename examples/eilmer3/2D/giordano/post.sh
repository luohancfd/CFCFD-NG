#!/bin/bash

e3post.py --job=inf_cyl --slice-list="0,:,0,0;2,:,0,0" --output=stag-prof-50Pa-Blackman.data
gnuplot plot-prof.gplot
