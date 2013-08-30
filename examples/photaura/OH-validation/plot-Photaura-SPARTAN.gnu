set term postscript eps enhanced color "Helvetica,24"
set output "emission-Photaura-SPARTAN.eps"
set size 1.5,1.0
set xlabel "Wavelength, {/Symbol l} (nm)"
set ylabel "Emission coefficient, j_{/Symbol l} (W/m^3-sr-m)"
set xrange [305:320]
set autoscale y
# set logscale y
set grid
set key top right
# NOTE: converting spartan data file from W/m3-sr-cm-1 to W/m3-sr-m via j_lambda = j_eta * eta / lambda
plot 'coefficient_spectra_with_AF.txt' using 1:2 with lines lw 1 lt 1 t 'Photaura', \
     'SPARTAN-spectra/nu_IE_IA_69.txt' using (1/$1*1.0e7):($2*$1*$1*1.0e2) with lines lw 1 lt 1 lc 3 t 'Spartan'

set output "absorption-Photaura-SPARTAN.eps"
set ylabel "Absorption coefficient, {/Symbol k}_{/Symbol l} (1/m)"
set xrange [305:320]
set grid
set key top right
plot 'coefficient_spectra_with_AF.txt' using 1:4 with lines lw 1 lt 1 t 'Photaura', \
     'SPARTAN-spectra/nu_IE_IA_69.txt' using (1/$1*1.0e7):($3) with lines lw 1 lt 1 lc 3 t 'Spartan'

set output "intensity-Photaura-SPARTAN.eps"
set ylabel "Spectral intensity, I_{/Symbol l} (W/m^2-sr-m)"
set xrange [305:320]
set grid
set key top right
plot 'intensity_spectra.txt' using 1:2 with lines lw 1 lt 1 t 'Photaura', \
     'SPARTAN-spectra/intensity_spectra.txt' using 1:2 with lines lw 1 lt 1 lc 3 t 'Spartan'
