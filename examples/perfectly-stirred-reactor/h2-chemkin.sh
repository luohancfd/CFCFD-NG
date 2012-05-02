#!/bin/bash
# h2-chemkin.sh
# Exercise the perfectly-stirred-reactor simulation 
# for the h2-chemkin reaction scheme.
#
# PJ, 06-Sep-2011 adapted from Brendan's work.
# 
# Note that this script puts Brendan's species files (containing chemkin
# thermo data) over the top of the installed species files.  This is done
# so that the reactor simulations match the Chemkin data closely, 
# otherwise, there is a delay in the reactions.
# At the end, it restores the original species files.

echo "Adjust thermo data data files."
SPECIES=${HOME}/e3bin/species
if [ ! -e saved-species-files ]; then
    mkdir saved-species-files
fi
rsync -a ${SPECIES}/*.lua ./saved-species-files
echo "Have saved original species files."
rsync -av ${SPECIES}/species-with-chemkin-thermo-data/*.lua ${SPECIES}/
echo "Have put Chemkin thermo data in place."

INPUTFILES=${HOME}/cfcfd3/lib/gas/reaction-schemes/methane-combustion
cp ${INPUTFILES}/h2-chemkin.lua .
cp ${INPUTFILES}/h2-chemkin.inp .
gasfile h2-chemkin.inp thermally-perfect-h2-chemkin.lua
echo "Done with gas model setup, now run reactor simulations."

time ./reactor.py \
    --type=cvfm \
    --gas=thermally-perfect-h2-chemkin.lua \
    --chem=h2-chemkin.lua \
    --T=1000 \
    --p=1 \
    --X="{'H2':1.0, 'O2':3.0, 'N2':0.1}" \
    --tmax=3.0e-04 --dt=1e-7 > cvfm-output.dat

time ./reactor.py \
    --type=cpfm \
    --gas=thermally-perfect-h2-chemkin.lua \
    --chem=h2-chemkin.lua \
    --T=1000 \
    --p=1 \
    --X="{'H2':1.0, 'O2':3.0, 'N2':0.1}" \
    --tmax=3e-04 --dt=1e-7 > cpfm-output.dat

time ./simple_reactor.py \
    --type=cvfm \
    --gas=thermally-perfect-h2-chemkin.lua \
    --chem=h2-chemkin.lua \
    --T=1000 \
    --p=1 \
    --X="{'H2':1.0, 'O2':3.0, 'N2':0.1}" \
    --tmax=3.0e-04 --dt=1e-7 > cvfm-sp-output.dat

time ./simple_reactor.py \
    --type=cpfm \
    --gas=thermally-perfect-h2-chemkin.lua \
    --chem=h2-chemkin.lua \
    --T=1000 \
    --p=1 \
    --X="{'H2':1.0, 'O2':3.0, 'N2':0.1}" \
    --tmax=3e-04 --dt=1e-7 > cpfm-sp-output.dat

gnuplot <<EOF
set terminal postscript eps 20
set output "cpfm-h2-chemkin-1000K-1bar-T.eps"
set title "h2-chemkin"
#
set mxtics 5
set xtics 0,1e-4,3e-4
#set xrange[0:1e-2]
set format x '%1.0e'
set xlabel 't, s'
#
set mytics 5
set ytics 1000,200,2800
set yrange[1000:2600]
set ylabel 'T, K'
#
set key bottom right
#
plot \
'h2-chemkin-chemkin-output-1000K-1bar.dat' u 1:2 w l t 'Chemkin-II', \
'cpfm-sp-output.dat' u 1:4 w l t 'simple_reactor', \
'cpfm-output.dat' u 1:5 w l t 'o/s simple_reactor'
EOF

gnuplot <<EOF
set terminal postscript eps 20
set output "cpfm-h2-chemkin-1000K-1bar-X.eps"
set title "h2-chemkin"
#
set mxtics 5
set xtics 0,1e-4,3e-4
set format x '%1.0e'
set xlabel 't, s'
#
set logscale y
set mytics 5 
set ytics mirror
set yrange[1e-12:1]
set ylabel 'X'
#
set key bottom right
#
# t, T, p, H2, H, O2, O, OH, HO2, H2O2, H2O, N, N2, NO
plot 'h2-chemkin-chemkin-output-1000K-1bar.dat' u 1:3 w l lt 1 t 'Chemkin-II',\
     '' u 1:4 w l lt 1 notitle, \
     '' u 1:5 w l lt 1 notitle, \
     '' u 1:6 w l lt 1 notitle, \
     '' u 1:7 w l lt 1 notitle, \
     '' u 1:8 w l lt 1 notitle, \
     '' u 1:9 w l lt 1 notitle, \
     '' u 1:10 w l lt 1 notitle, \
     '' u 1:11 w l lt 1 notitle, \
     '' u 1:12 w l lt 1 notitle, \
     '' u 1:13 w l lt 1 notitle, \
     'cpfm-output.dat' u 1:3 w l lt 2 t 'simple_reactor', \
     '' u 1:4 w l lt 2 notitle, \
     '' u 1:5 w l lt 2 notitle, \
     '' u 1:6 w l lt 2 notitle, \
     '' u 1:7 w l lt 2 notitle, \
     '' u 1:8 w l lt 2 notitle, \
     '' u 1:9 w l lt 2 notitle, \
     '' u 1:10 w l lt 2 notitle, \
     '' u 1:11 w l lt 2 notitle, \
     '' u 1:12 w l lt 2 notitle, \
     'cpfm-sp-output.dat' u 1:7 w l lt 3 t 'o/s simple_reactor',\
     '' u 1:8 w l lt 3 notitle, \
     '' u 1:9 w l lt 3 notitle, \
     '' u 1:10 w l lt 3 notitle, \
     '' u 1:11 w l lt 3 notitle, \
     '' u 1:12 w l lt 3 notitle, \
     '' u 1:13 w l lt 3 notitle, \
     '' u 1:14 w l lt 3 notitle, \
     '' u 1:15 w l lt 3 notitle, \
     '' u 1:16 w l lt 3 notitle, \
     '' u 1:17 w l lt 3 notitle
EOF

echo "Restore the original species files."
rsync -av ./saved-species-files/*.lua ${SPECIES}/
echo "Have restored original species files."
