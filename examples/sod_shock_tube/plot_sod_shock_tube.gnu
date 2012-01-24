set terminal epslatex
set output 'sod_profile.tex'

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
set yrange[0:1.2e5]
set format y '%2.1t\e{%L}'
set ylabel '$p$, Pa'

set key top right

plot \
'sod-analytical.dat' u 1:3 w l lt 1 t 'ideal gas proceedure', \
'sod-numerical.dat' u 1:3 w l lt 2 t 'real gas proceedure'
