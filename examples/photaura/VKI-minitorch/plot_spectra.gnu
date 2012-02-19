set term postscript eps enhanced color "Helvetica,24"
set output "Intensity3500-4300A.eps"
set size 1.5,1.0
set xlabel "Wavelength, {/Symbol l} (Ang)"
set ylabel "Spectral intensity, I\_{{/Symbol l}} (W/m^2-sr-{/Symbol m}m)"
set xrange [3500:4300]
set autoscale y
# set logscale y
set grid
set key top right
plot 'plasma-transported-intensity.txt' using ($1*10):($2*1.0e-6) with lines lw 1 lt 1 t 'Calculated', \
     'plasma-transported-intensity-optically-thin.txt' using ($1*10):($2*1.0e-6) with lines lw 1 lt 3 t 'Calculated (optically thin)', \
     'exp_data/Intensity3500-4300A.dat' using ($1):($2) with lines lw 2 lt -1 t 'Measured'
     
set term postscript eps enhanced color "Helvetica,24"
set output "Intensity4300-10000A.eps"
set size 1.5,1.0
set xlabel "Wavelength, {/Symbol l} (Ang)"
set ylabel "Spectral intensity, I\_{{/Symbol l}} (W/m^2-sr-{/Symbol m}m)"
set xrange [4300:10000]
set autoscale y
# set logscale y
set grid
set key top left
plot 'plasma-transported-intensity.txt' using ($1*10):($2*1.0e-6) with lines lw 1 lt 1 t 'Calculated', \
     'plasma-transported-intensity-optically-thin.txt' using ($1*10):($2*1.0e-6) with lines lw 1 lt 3 t 'Calculated (optically thin)', \
     'exp_data/Intensity4300-10000A.dat' using ($1):($2) with lines lw 2 lt -1 t 'Measured'
     
set term postscript eps enhanced color "Helvetica,24"
set output "Emissivity-r1.24mm.eps"
set size 1.5,1.0
set xlabel "Wavelength, {/Symbol l} (Ang)"
set ylabel "Spectral emissivity, I\_{{/Symbol l}} (W/m^2-sr-{/Symbol m}m)"
set xrange [4300:10000]
set autoscale y
# set logscale y
set grid
set key top left
plot 'emissivities_at_1240_microns.txt' using ($1*10):($2*1.0e-6*0.71) with lines lw 1 lt 1 t 'Calculated x 0.71', \
     'exp_data/Emissivity-r1.24mm.dat' using ($1):($2) with lines lw 2 lt -1 t 'Measured'
     
set term postscript eps enhanced color "Helvetica,24"
set output "Emissivity-r3.82mm.eps"
set size 1.5,1.0
set xlabel "Wavelength, {/Symbol l} (Ang)"
set ylabel "Spectral emissivity, I\_{{/Symbol l}} (W/m^2-sr-{/Symbol m}m)"
set xrange [4300:9000]
set autoscale y
# set logscale y
set grid
set key top left
plot 'emissivities_at_3820_microns.txt' using ($1*10):($2*1.0e-6*0.62) with lines lw 1 lt 1 t 'Calculated x 0.62', \
     'exp_data/Emissivity-r3.82mm.dat' using ($1):($2) with lines lw 2 lt -1 t 'Measured'
     
set term postscript eps enhanced color "Helvetica,24"
set output "Emissivity-r8.43mm.eps"
set size 1.5,1.0
set xlabel "Wavelength, {/Symbol l} (Ang)"
set ylabel "Spectral emissivity, I\_{{/Symbol l}} (W/m^2-sr-{/Symbol m}m)"
set xrange [4300:9000]
set autoscale y
# set logscale y
set grid
set key top left
plot 'emissivities_at_8430_microns.txt' using ($1*10):($2*1.0e-6*0.5) with lines lw 1 lt 1 t 'Calculated x 0.5', \
     'exp_data/Emissivity-r8.43mm.dat' using ($1):($2) with lines lw 2 lt -1 t 'Measured'
