# plot_data.gnu

set term postscript eps enhanced 20
set output "profile_moles.eps"
set title "Species profiles behind normal shock - Titan gas"
set xlabel "x, m"
set ylabel "mole fraction"
set logscale y
set key bottom left
plot 'ballute.data' u 1:9 t 'N_2' w l , \
     'ballute.data' u 1:10 t 'CH4' w l , \
     'ballute.data' u 1:11 t 'CH3' w l , \
     'ballute.data' u 1:15 t 'H2' w l, \
     'ballute.data' u 1:16 t 'CN' w l, \
     'ballute.data' u 1:19 t 'N' w l, \
     'ballute.data' u 1:20 t 'C' w l
