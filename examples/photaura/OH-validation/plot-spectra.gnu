set term postscript eps enhanced color "Helvetica,24"
set output "emission-spectra.eps"
set size 1.5,1.0
set xlabel "Wavelength, nm"
set ylabel "Emission coefficient, j_{/Symbol l} (arbitrary units)"
set xrange [305:320]
set autoscale y
# set logscale y
set grid
set key top right
# Use STATS function in gnuplot to find maximum in Photaura's intensity 
# spectra for normalisation to arbitrary units.
stats 'coefficient_spectra_with_AF.txt' using 2
photaura_max=STATS_max
stats 'SPARTAN-spectra/nu_IE_IA_69.txt' using ($2*$1*$1/1.0e7)
spartan_max=STATS_max
plot 'coefficient_spectra_with_AF.txt' using ($1):(($2)*100/photaura_max) with lines lw 1 lt 1 t 'Photaura', \
     'LIFBASE-spectra/OH_AX_res_2_5_3000K.mod' using ($1/10):($2) with lines lw 1 lt 1 lc 3 t 'LIFBASE', \
     'SPARTAN-spectra/nu_IE_IA_69.txt' using (1/$1*1.0e7):($2*$1*$1/1.0e7/spartan_max*100) with lines lw 1 lt 1 lc 4 t 'SPARTAN'

