set term postscript eps enhanced color "Helvetica,22"
set grid
set size 1.25,1.0

set output "intensity_spectra.eps"
set xlabel "Wavelength, {/Symbol l} (nm)"
set ylabel "Spectral radiance, I_{/Symbol l} (W/cm^2-sr-{/Symbol m}m)"
set logscale y
set key bottom
#set yrange [1:1e5]
plot      'NIST_intensity_spectra.txt' u 1:($2*1.0e-10) w l lt -1 lw 1 t 'NIST', \
          'TOPBase_intensity_spectra.txt' u 1:($2*1.0e-10) w l lt 1 lw 1 t 'TOPBase'

set output "intensity_spectra_with_cumulative.eps"
set xlabel "Wavelength, {/Symbol l} (nm)"
set ylabel "Spectral radiance, I_{/Symbol l} (W/cm^2-sr-{/Symbol m}m)"
set logscale y
set y2label "Cumulative radiance, I (W/cm^2-sr)"
set key bottom
set y2tics
#set yrange [1:1e5]
plot      'NIST_intensity_spectra.txt' u 1:($2*1.0e-10) w l lt -1 lw 1 t 'NIST', \
          'NIST_intensity_spectra.txt' u 1:($3*1.0e-4) axes x1y2 w l lt -1 lw 2 notitle, \
          'TOPBASE_intensity_spectra.txt' u 1:($2*1.0e-10) w l lt 1 lw 1 t 'TOPBase', \
	  'TOPBASE_intensity_spectra.txt' u 1:($3*1.0e-4) axes x1y2 w l lt 1 lw 2 notitle
          
set autoscale y
set autoscale y2
          
set output "emission_spectra_with_cumulative.eps"
set xlabel "Wavelength, {/Symbol l} (nm)"
set ylabel "Spectral emission, j_{/Symbol l} (W/cm^3-sr-{/Symbol m}m)"
set logscale y
set y2label "Cumulative emission, j (W/cm^3-sr)"
set key bottom
set y2tics 
#set yrange [0.01:1e8]
plot      'NIST_coefficient_spectra.txt' u 1:($2*1.0e-12) w l lt -1 lw 1 t 'NIST', \
          'NIST_coefficient_spectra.txt' u 1:($3*1.0e-6) axes x1y2 w l lt -1 lw 3 notitle, \
	  'TOPBase_coefficient_spectra.txt' u 1:($2*1.0e-12) w l lt 1 lw 1 t 'TOPBase', \
          'TOPBase_coefficient_spectra.txt' u 1:($3*1.0e-6) axes x1y2 w l lt 1 lw 3 notitle
