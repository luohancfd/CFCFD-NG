#!/bin/bash
# drm19.sh
# Exercise the perfectly-stirred-reactor simulation 
# for the DRM19 reaction scheme.

INPUTFILES=${HOME}/cfcfd2/lib/gas/reaction-schemes/methane-combustion
cp ${INPUTFILES}/drm19.lua .
cp ${INPUTFILES}/drm19.inp .
gasfile drm19.inp thermally-perfect-drm19.lua

time ./simple_reactor.py \
    --type=cvfm \
    --gas=thermally-perfect-drm19.lua \
    --chem=drm19.lua \
    --T=2000 \
    --p=1 \
    --X="{'CH4':1.0, 'O2':2.0, 'N2':7.52}" \
    --tmax=4.0e-04 > cvfm-sp-output.dat

time ./simple_reactor.py \
    --type=cpfm \
    --gas=thermally-perfect-drm19.lua \
    --chem=drm19.lua \
    --T=2000 \
    --p=1 \
    --X="{'CH4':1.0, 'O2':2.0, 'N2':7.52}" \
    --tmax=4e-04 > cpfm-sp-output.dat

gnuplot <<EOF
set terminal postscript eps 20
set output "cvfm-drm19-2000K-1bar-T.eps"
#
set mxtics 5
set xtics 0,1e-4,4e-4
#set xrange[0:1e-2]
set format x '%1.0e'
set xlabel 't, s'
#
set mytics 5
set ytics 2000,200,3200
set yrange[2000:3200]
set ylabel 'T, K'
#
set key bottom right
#
plot 'cvfm-sp-output.dat' u 1:4 w l lt 1 t 'simple_reactor w/ DRM19'
EOF

gnuplot <<EOF
set terminal postscript eps 20
set output "cvfm-drm19-2000K-1bar-X.eps"
#
set mxtics 5
set xtics 0,1e-4,4e-4
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
plot 'cvfm-sp-output.dat' u 1:7 w l lt 1 t 'simple_reactor w/ DRM19', \
    '' u 1:8 w l lt 1 notitle, \
    '' u 1:9 w l lt 1 notitle, \
    '' u 1:10 w l lt 1 notitle, \
    '' u 1:11 w l lt 1 notitle, \
    '' u 1:12 w l lt 1 notitle, \
    '' u 1:13 w l lt 1 notitle, \
    '' u 1:14 w l lt 1 notitle, \
    '' u 1:15 w l lt 1 notitle, \
    '' u 1:16 w l lt 1 notitle, \
    '' u 1:17 w l lt 1 notitle, \
    '' u 1:18 w l lt 1 notitle, \
    '' u 1:19 w l lt 1 notitle, \
    '' u 1:20 w l lt 1 notitle, \
    '' u 1:21 w l lt 1 notitle, \
    '' u 1:22 w l lt 1 notitle, \
    '' u 1:23 w l lt 1 notitle, \
    '' u 1:24 w l lt 1 notitle, \
    '' u 1:25 w l lt 1 notitle, \
    '' u 1:26 w l lt 1 notitle, \
    '' u 1:27 w l lt 1 notitle
EOF

gnuplot <<EOF
set terminal postscript eps 20
set output "cpfm-drm19-2000K-1bar-T.eps"
#
set mxtics 5
set xtics 0,1e-4,4e-4
#set xrange[0:1e-2]
set format x '%1.0e'
set xlabel 't, s'
#
set mytics 5
set ytics 2000,200,3200
set yrange[2000:3200]
set ylabel 'T, K'

set key bottom right

plot 'drm19-chemkin-output-2000K-1bar.dat' u 1:2 w l lt 1 t 'Chemkin-II w/ DRM19', \
     'cpfm-sp-output.dat' u 1:4 w l lt 2 t 'o/s simple_reactor w/ DRM19'
EOF

gnuplot <<EOF
set terminal postscript eps 20
set output "cpfm-drm19-2000K-1bar-X.eps"
#
set mxtics 5
set xtics 0,1e-4,4e-4
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
plot 'drm19-chemkin-output-2000K-1bar.dat' u 1:3 w l lt 1 t 'Chemkin-II', \
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
    '' u 1:14 w l lt 1 notitle, \
    '' u 1:15 w l lt 1 notitle, \
    '' u 1:16 w l lt 1 notitle, \
    '' u 1:17 w l lt 1 notitle, \
    '' u 1:18 w l lt 1 notitle, \
    '' u 1:19 w l lt 1 notitle, \
    '' u 1:20 w l lt 1 notitle, \
    '' u 1:21 w l lt 1 notitle, \
    '' u 1:22 w l lt 1 notitle, \
    '' u 1:23 w l lt 1 notitle, \
    'cpfm-sp-output.dat' u 1:7 w l lt 2 t 'simple_reactor w/ DRM19', \
    '' u 1:8 w l lt 2 notitle, \
    '' u 1:9 w l lt 2 notitle, \
    '' u 1:10 w l lt 2 notitle, \
    '' u 1:11 w l lt 2 notitle, \
    '' u 1:12 w l lt 2 notitle, \
    '' u 1:13 w l lt 2 notitle, \
    '' u 1:14 w l lt 2 notitle, \
    '' u 1:15 w l lt 2 notitle, \
    '' u 1:16 w l lt 2 notitle, \
    '' u 1:17 w l lt 2 notitle, \
    '' u 1:18 w l lt 2 notitle, \
    '' u 1:19 w l lt 2 notitle, \
    '' u 1:20 w l lt 2 notitle, \
    '' u 1:21 w l lt 2 notitle, \
    '' u 1:22 w l lt 2 notitle, \
    '' u 1:23 w l lt 2 notitle, \
    '' u 1:24 w l lt 2 notitle, \
    '' u 1:25 w l lt 2 notitle, \
    '' u 1:26 w l lt 2 notitle, \
    '' u 1:27 w l lt 2 notitle
EOF
