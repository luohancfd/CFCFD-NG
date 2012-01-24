

set term postscript eps enhanced color 24
set size 2.1,1.5
set logscale y
set output "140-900nm-intensity.eps"
set xlabel "Wavelength, {/Symbol l} (nm)" font "Helvetica,24"
set ylabel "Emission intensity, I\_{{/Symbol l}} (W/cm^2-{/Symbol m}m-sr)" font "Helvetica,24"
set grid
set key top right
plot '2T-basic-Lee-intensity-at-peak.dat' using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  3 t '2T Basic Lee', \
     '2T-full-PLF-intensity-at-peak.dat'      using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  1 t '2T Full PLF', \
     '2T-full-PLFM-intensity-at-peak.dat'     using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  2 t '2T Full PLFM', \
     '3T-full-PLF-intensity-at-peak.dat'      using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  4 t '3T Full PLF', \
     'E090806_IvW.txt'                    using 1:2 every ::228::820 w l lw 2 lt -1 t 'Combined EAST Shot Data', \
     'E091206_IvW.txt'                    using 1:2 every ::167::820 w l lw 2 lt -1 notitle, \
     'W090806_IvW.txt'                    using 1:2 w l lw 2 lt -1 notitle, \
     'W091206_IvW.txt'                    using 1:2 w l lw 2 lt -1 notitle
     
set output "440-580nm-intensity.eps"
set xlabel "Wavelength, {/Symbol l} (nm)" font "Helvetica,24"
set ylabel "Emission intensity, I\_{{/Symbol l}} (W/cm^2-{/Symbol m}m-sr)" font "Helvetica,24"
set xrange [440:580]
set yrange [1.0e-2:1.0e3]
set grid
set key top right
plot '2T-basic-Lee-intensity-at-peak.dat' using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  3 t '2T Basic Lee', \
     '2T-full-PLF-intensity-at-peak.dat'  using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  1 t '2T Full PLF', \
     '2T-full-PLFM-intensity-at-peak.dat' using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  2 t '2T Full PLFM', \
     'W091206_IvW.txt'                    using 1:2 w l lw 2 lt -1 t 'EAST Shot 46/19 IR', \
     'E091206_IvW.txt'                    using 1:2 every ::167::820 w l lw 1 lt -1 t 'EAST Shot 46/19 UV'
     
set term postscript eps enhanced color 24
set size 2.1,1.5
set output "340-580nm-intensity-lin.eps"
set xlabel "Wavelength, {/Symbol l} (nm)" font "Helvetica,24"
set ylabel "Emission intensity, I\_{{/Symbol l}} (W/cm^2-{/Symbol m}m-sr)" font "Helvetica,24"
set xrange [340:580]
set yrange [0:1200]
set grid
set key top right
plot '2T-full-Lee-intensity-at-peak.dat'   using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  3 t '2T Full Lee', \
     '2T-full-Park-intensity-at-peak.dat'  using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  1 t '2T Full Park', \
     '3T-full-PLF-intensity-at-peak.dat'   using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  4 t '3T Full PLF', \
     'W091206_IvW.txt'                    using 1:2 w l lw 2 lt -1 t 'EAST Shot 46/19 Combined IR and UV', \
     'E091206_IvW.txt'                    using 1:2 w l lw 2 lt -1 notitle

set term postscript eps enhanced color 28
set size 2.1,1.5
set logscale y
set output "180-320nm-intensity.eps"
set xlabel "Wavelength, {/Symbol l} (nm)" font "Helvetica,24"
set ylabel "Emission intensity, I\_{{/Symbol l}} (W/cm^2-{/Symbol m}m-sr)" font "Helvetica,24"
set xrange [180:320]
set yrange [1.0e-1:1.0e4]
set grid
set key top right
plot '2T-basic-Lee-intensity-at-peak.dat' using 1:($2*1.0e-4*1.0e-6) with lines lw 6 lt  3 t '2T Basic Lee', \
     '2T-full-Park-intensity-at-peak.dat' using 1:($2*1.0e-4*1.0e-6) with lines lw 6 lt  4 t '2T Full Park', \
     '2T-full-PLF-intensity-at-peak.dat'      using 1:($2*1.0e-4*1.0e-6) with lines lw 6 lt  1 t '2T Full PLF', \
     '2T-full-PLFM-intensity-at-peak.dat'     using 1:($2*1.0e-4*1.0e-6) with lines lw 6 lt  2 t '2T Full PLFM', \
     'E090806_IvW.txt'                    using 1:2 every ::230 w l lw 3 lt -1 t 'EAST Shot 46/18 UV'
     
set output "700-900nm-intensity-lin.eps"
set xlabel "Wavelength, {/Symbol l} (nm)" font "Helvetica,24"
set ylabel "Emission intensity, I\_{{/Symbol l}} (W/cm^2-{/Symbol m}m-sr)" font "Helvetica,24"
set xrange [700:900]
set yrange [0.1:2000]
set grid
set key top right
plot '2T-basic-Lee-intensity-at-peak.dat' using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  3 t '2T Basic Lee', \
     '2T-full-PLF-intensity-at-peak.dat'      using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  1 t '2T Full PLF', \
     '2T-full-PLFM-intensity-at-peak.dat'     using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  2 t '2T Full PLFM', \
     '3T-full-PLF-intensity-at-peak.dat'      using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  4 t '3T Full PLF', \
     'W090806_IvW.txt'                    using 1:2 w l lw 2 lt -1 t 'EAST Shot 46/18 UV'
     
unset logscale y

set output "440-580nm-intensity-lin.eps"
set xlabel "Wavelength, {/Symbol l} (nm)" font "Helvetica,24"
set ylabel "Emission intensity, I\_{{/Symbol l}} (W/cm^2-{/Symbol m}m-sr)" font "Helvetica,24"
set xrange [440:580]
set yrange [0:500]
set grid
set key top right
plot '2T-full-Lee-intensity-at-peak.dat'   using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  3 t '2T Full Lee', \
     '2T-full-PLF-intensity-at-peak.dat'   using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  1 t '2T Full PLF', \
     '2T-full-PLFM-intensity-at-peak.dat'  using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  2 t '2T Full PLFM', \
     '3T-full-PLF-intensity-at-peak.dat'   using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  4 t '3T Full PLF', \
     'W091206_IvW.txt'                    using 1:2 w l lw 2 lt -1 t 'EAST Shot 46/19 IR', \
     'E091206_IvW.txt'                    using 1:2 w l lw 1 lt -1 t 'EAST Shot 46/19 UV'
     
set term postscript eps enhanced color 30
set size 2.1,1.5
set output "340-580nm-intensity.eps"
set xlabel "Wavelength, {/Symbol l} (nm)" font "Helvetica,30"
set ylabel "Intensity, I\_{{/Symbol l}} (W/cm^2-{/Symbol m}m-sr)" font "Helvetica,30"
set xrange [340:580]
set yrange [1e0:1200]
set logscale y
set grid
set key top right
plot '2T-full-Lee-intensity-at-peak.dat'   using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  3 t '2T Full Lee', \
     '2T-full-Park-intensity-at-peak.dat'  using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  1 t '2T Full Park', \
     '3T-full-PLF-intensity-at-peak.dat'   using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  4 t '3T Full PLF', \
     'W091206_IvW.txt'                    using 1:2 w l lw 2 lt -1 t 'EAST Shot 46/19 Combined IR and UV', \
     'E091206_IvW.txt'                    using 1:2 every ::::840 w l lw 2 lt -1 notitle
     
set output "700-900nm-intensity-lin.eps"
set xlabel "Wavelength, {/Symbol l} (nm)" font "Helvetica,24"
set ylabel "Emission intensity, I\_{{/Symbol l}} (W/cm^2-{/Symbol m}m-sr)" font "Helvetica,24"
set xrange [700:900]
set yrange [0.0:1000]
set grid
set key top right
plot '2T-basic-Lee-intensity-at-peak.dat' using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  3 t '2T Basic Lee', \
     '2T-full-PLF-intensity-at-peak.dat'      using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  1 t '2T Full PLF', \
     '2T-full-PLFM-intensity-at-peak.dat'     using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  2 t '2T Full PLFM', \
     '3T-full-PLF-intensity-at-peak.dat'      using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  4 t '3T Full PLF', \
     'W090806_IvW.txt'                    using 1:2 w l lw 2 lt -1 t 'EAST Shot 46/18 UV'
