#!/bin/bash
# run-and-plot.sh
# Exercise the perfectly-stirred-reactor simulation 
# for the GRI-Mech3.0 reaction scheme.
#
# PJ, 05-Sep-2011

time e3shared.exe --job=psr --run
e3history.py --ijk="0,0,0" --output="eilmer3-history.data" \
    --add-molef --gmodel-file="thermally-perfect-grimech30.lua" \
    hist/psr.hist.b0000

gnuplot <<EOF
set terminal postscript eps 20
set output "cvfm-grimech30-2000K-1bar-T.eps"
#
set title "Methane combustion with GRI-Mech3.0 reaction scheme."
set mxtics 5
set xtics 0,100,400
set xlabel 't, microseconds'
#
set mytics 5
set ytics 2000,200,3200
set yrange[2000:3200]
set ylabel 'T, K'
#
set key bottom right
#
plot 'cvfm-sp-output.dat' using (\$1*1.0e6):4 with lines linetype 1 title 'simple_reactor', \
     'eilmer3-history.data' using (\$1*1.0e6):77 with lines linetype 2 title 'Eilmer3'
EOF


gnuplot <<EOF
set terminal postscript eps 20
set output "cvfm-grimech30-2000K-1bar-p.eps"
#
set title "Methane combustion with GRI-Mech3.0 reaction scheme."
set mxtics 5
set xtics 0,100,400
set xlabel 't, microseconds'
#
set mytics 5
set ytics 1,0.2,2
set yrange[1:2]
set ylabel 'p, bar'
#
set key bottom right
#
plot 'cvfm-sp-output.dat' using (\$1*1.0e6):(\$5/1.0e5) with lines linetype 1 title 'simple_reactor', \
     'eilmer3-history.data' using (\$1*1.0e6):(\$13/1.0e5) with lines linetype 2 title 'Eilmer3'
EOF

gnuplot <<EOF
set terminal postscript eps 20
set output "cvfm-grimech30-2000K-1bar-X-H2O.eps"
#
set title "Methane combustion with GRI-Mech3.0 reaction scheme."
set mxtics 5
set xtics 0,100,400
set xlabel 't, microseconds'
#
set logscale y
set mytics 5 
set ytics mirror
set yrange[1e-12:1]
set ylabel 'mole fraction H2O'
#
set key bottom right
#
plot 'cvfm-sp-output.dat' using (\$1*1.0e6):12 with lines linetype 1 title 'simple_reactor', \
     'eilmer3-history.data' using (\$1*1.0e6):83 with lines linetype 2 title 'Eilmer3'
EOF
