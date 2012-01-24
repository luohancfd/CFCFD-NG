set terminal epslatex
set output 'real_sod_profile-u.tex'

# set x-axis
#set logscale x
set mxtics 5
set xtics
#set xrange[1e3:1e9]
set format x '%2.1f'
set xlabel '$x$, m'

# set y-axis
#set logscale y
set mytics 5
set yrange[0:250]
set format y '%2.1t\e{%L}'
set ylabel '$u$, m.s\pow{-1}'

set key top left

plot \
'sod-tpg.dat' u 1:5 w l lt 1 t 'thermally perfect', \
'sod-nag.dat' u 1:5 w l lt 2 t 'Noble-Abel',\
'sod-vdw.dat' u 1:5 w l lt 3 t 'van der Waals'
