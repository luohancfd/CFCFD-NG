set term postscript eps enhanced color "Helvetica,22"
set grid
set size 1.25,1.0

set output "intensity_spectra.eps"
set xlabel "Wavelength, {/Symbol l} (nm)"
set ylabel "Spectral radiance, I_{/Symbol l} (W/cm^2-sr-{/Symbol m}m)"
set logscale y
set key bottom
set yrange [1:1e5]
plot      'legacy_intensity_spectra.txt' u 1:($2*1.0e-10) w l lt -1 lw 1 t 'Legacy spectral modelling', \
          'new_intensity_spectra.txt' u 1:($2*1.0e-10) w l lt 1 lw 1 t 'Advanced spectral modelling'

set output "intensity_spectra_with_cumulative.eps"
set xlabel "Wavelength, {/Symbol l} (nm)"
set ylabel "Spectral radiance, I_{/Symbol l} (W/cm^2-sr-{/Symbol m}m)"
set logscale y
set y2label "Cumulative radiance, I (W/cm^2-sr)"
set key bottom
set y2tics 20
set yrange [1:1e5]
plot      'legacy_intensity_spectra.txt' u 1:($2*1.0e-10) w l lt -1 lw 1 t 'Legacy spectral modelling', \
          'legacy_intensity_spectra.txt' u 1:($3*1.0e-4) axes x1y2 w l lt -1 lw 2 notitle, \
          'new_intensity_spectra.txt' u 1:($2*1.0e-10) w l lt 1 lw 1 t 'Advanced spectral modelling', \
          'new_intensity_spectra.txt' u 1:($3*1.0e-4) axes x1y2 w l lt 1 lw 2 notitle