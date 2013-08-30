set term postscript eps enhanced color "Helvetica,24"
set xlabel "Wavelength, {/Symbol l} (nm)"
set ylabel "Spectral intensity, I_{/Symbol l} (arbitrary units)"
set xrange [305:320]
set autoscale y
set grid

set output "intensity-Photaura-exp-SL1.eps"
set size 1.5,1.0
set key top right
# Use STATS function in gnuplot to find maximum in Photaura's intensity 
# spectra for normalisation to arbitrary units.
stats 'EXP-spectra/experimental_data1.txt' using 2
SL1_I_max = STATS_max
stats 'intensity_spectra.txt' using 2
I_max = STATS_max
stats 'coefficient_spectra_with_AF.txt' using 2
j_max = STATS_max
plot 'EXP-spectra/experimental_data1.txt' using  ($1):($2*100/SL1_I_max) with lines lw 1 lt 1 lc 2 t 'EXP-slit-loc1', \
     'coefficient_spectra_with_AF.txt' using ($1):(($2)*100/j_max) with lines lw 1 lt 1 lc 3 t 'Photaura: Optically thin'

set output "intensity-Photaura-exp-SL2.eps"
set key top right
# Use STATS function in gnuplot to find maximum in Photaura's intensity 
# spectra for normalisation to arbitrary units.
stats 'EXP-spectra/experimental_data2.txt' using 2
SL2_I_max = STATS_max
stats 'intensity_spectra.txt' using 2
I_max = STATS_max
stats 'coefficient_spectra_with_AF.txt' using 2
j_max = STATS_max
plot 'EXP-spectra/experimental_data2.txt' using  ($1):($2*100/SL2_I_max) with lines lw 1 lt 1 lc 4 t 'EXP-slit-loc2', \
     'coefficient_spectra_with_AF.txt' using ($1):(($2)*100/j_max) with lines lw 1 lt 1 lc 3 t 'Photaura: Optically thin'

set output "intensity-Photaura-exp-SL3.eps"
set key bottom right
# Use STATS function in gnuplot to find maximum in Photaura's intensity 
# spectra for normalisation to arbitrary units.
stats 'EXP-spectra/experimental_data3.txt' using 2
SL3_I_max = STATS_max
stats 'intensity_spectra.txt' using 2
I_max = STATS_max
stats 'coefficient_spectra_with_AF.txt' using 2
j_max = STATS_max
plot 'EXP-spectra/experimental_data3.txt' using  ($1):($2*100/SL3_I_max) with lines lw 1 lt 1 lc 7 t 'EXP-slit-loc3', \
     'intensity_spectra.txt' using ($1):(($2)*100/I_max) with lines lw 1 lt 1 lc 1 t 'Photaura: Absorption over {/Symbol D}L=10cm'

set output "intensity-Photaura-exp.eps"
set size 1.5,1.5
set key below
# Use STATS function in gnuplot to find maximum in Photaura's intensity 
# spectra for normalisation to arbitrary units.
stats 'EXP-spectra/experimental_data1.txt' using 2
SL1_I_max = STATS_max
stats 'EXP-spectra/experimental_data2.txt' using 2
SL2_I_max = STATS_max
stats 'EXP-spectra/experimental_data3.txt' using 2
SL3_I_max = STATS_max
stats 'intensity_spectra.txt' using 2
I_max = STATS_max
stats 'coefficient_spectra_with_AF.txt' using 2
j_max = STATS_max
plot 'EXP-spectra/experimental_data1.txt' using  ($1):($2*100/SL1_I_max) with lines lw 1 lt 1 lc 2 t 'EXP-slit-loc1', \
     'EXP-spectra/experimental_data2.txt' using  ($1):($2*100/SL2_I_max) with lines lw 1 lt 1 lc 4 t 'EXP-slit-loc2', \
     'EXP-spectra/experimental_data3.txt' using  ($1):($2*100/SL3_I_max) with lines lw 1 lt 1 lc 7 t 'EXP-slit-loc3', \
     'intensity_spectra.txt' using ($1):(($2)*100/I_max) with lines lw 1 lt 1 lc 1 t 'Photaura: Absorption over {/Symbol D}L=10cm', \
     'coefficient_spectra_with_AF.txt' using ($1):(($2)*100/j_max) with lines lw 1 lt 1 lc 3 t 'Photaura: Optically thin'

