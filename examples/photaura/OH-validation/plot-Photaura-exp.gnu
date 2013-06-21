set term postscript eps enhanced color "Helvetica,24"
set output "Intensity-spectra-Photaura-exp.eps"
set size 1.5,1.0
set xlabel "Wavelength, nm"
set ylabel "Spectral intensity, I_{/Symbol l} (arbitrary units)"
set xrange [300:330]
set autoscale y
# set logscale y
set grid
set key top right
# Use STATS function in gnuplot to find maximum in Photaura's intensity 
# spectra for normalisation to arbitrary units.
stats 'intensity_spectra.txt' using 2
plot 'intensity_spectra.txt' using ($1):(($2)*100/STATS_max) with lines lw 1 lt 1 t 'Photaura', \
     'EXP-spectra/experimental_data1.txt' using  ($1):($2) with lines lw 1 lt 1 lc 2 t 'EXP-slit-loc1','EXP-spectra/experimental_data2.txt' using  ($1):($2) with lines lw 1 lt 1 lc 4 t 'EXP-slit-loc2','EXP-spectra/experimental_data3.txt' using  ($1):($2) with lines lw 1 lt 1 lc 7 t 'EXP-slit-loc3'


