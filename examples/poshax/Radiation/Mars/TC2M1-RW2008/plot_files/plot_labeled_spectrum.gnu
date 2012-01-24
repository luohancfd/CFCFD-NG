
set term postscript eps enhanced color 24
set size 2.1,1.5
set logscale y
set output "labeled-spectrum.eps"
set xlabel "Wavelength, {/Symbol l} (nm)" font "Helvetica,24"
set ylabel "Intensity, I\_{{/Symbol l}} (W/cm^2-{/Symbol m}m-sr)" font "Helvetica,24"
set grid
set key top right
set label "     NO {/Symbol g}, {/Symbol b}, {/Symbol d}, {/Symbol e}\n 215.2 +/- 24.8 nm" at 160, 4e4
set arrow from 215.2, 1.5e4 to 190.4, 1.5e4 lw 1
set arrow from 215.2, 1.5e4 to 240, 1.5e4 lw 1
set label "       C line\n 247.65 +/- 5 nm" at 200, 4e3
set arrow from 247.65, 1.5e3 to 242.65, 1.5e3 lw 2
set arrow from 247.65, 1.5e3 to 252.65, 1.5e3 lw 2
# set label "CN Violet {/Symbol D}v=0\n381 +/- 10 nm" at 310, 4e3
# set arrow from 381, 1.5e3 to 371, 1.5e3 lw 2
# set arrow from 381, 1.5e3 to 391, 1.5e3 lw 2
set label "CN Violet {/Symbol D}v=-1,0\n   366 +/- 25 nm" at 310, 4e4
set arrow from 366, 1.5e4 to 341, 1.5e4 lw 1
set arrow from 366, 1.5e4 to 391, 1.5e4 lw 1
set label "C_2 Swan {/Symbol D}v=-2,-1\n452.5 +/- 27.5 nm" at 390, 4.0e3
set arrow from 452.5, 1.5e3 to 425, 1.5e3 lw 1
set arrow from 452.5, 1.5e3 to 480, 1.5e3 lw 1
set label " C_2 Swan {/Symbol D}v=0,1\n527.5 +/- 42.5 nm" at 472.5, 2.0e4
set arrow from 527.5, 8.5e3 to 485, 8.5e3 lw 1
set arrow from 527.5, 8.5e3 to 570, 8.5e3 lw 1
# set label "      UV\n440 +/- 40 nm" at 400, 4.0e0
# set arrow from 440, 6e0 to 400, 6e0 lw 1
# set arrow from 440, 6e0 to 480, 6e0 lw 1
set label "   O triplet\n777.11 +/- 20 nm" at 727.11, 4e3
set arrow from 777.11, 1.5e3 to 767.11, 1.5e3 lw 1
set arrow from 777.11, 1.5e3 to 787.11, 1.5e3 lw 1
set yrange [1.0:1e5]
plot '2T-full-PLF-intensity-at-peak.dat'  using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  1 t '2T Full PLF computation', \
     'E090806_IvW.txt'                    using 1:2 every ::228::820 w l lw 2 lt -1 t 'Combined EAST Shot Data', \
     'E091206_IvW.txt'                    using 1:2 every ::167::840 w l lw 2 lt -1 notitle, \
     'W090806_IvW.txt'                    using 1:2 w l lw 2 lt -1 notitle, \
     'W091206_IvW.txt'                    using 1:2 w l lw 2 lt -1 notitle
     
unset label
unset arrow
     
set output "180-350nm-intensity.eps"
set xlabel "Wavelength, {/Symbol l} (nm)" font "Helvetica,24"
set ylabel "Emission intensity, I\_{{/Symbol l}} (W/cm^2-{/Symbol m}m-sr)" font "Helvetica,24"
set xrange [180:350]
set yrange [1.0e0:1.0e4]
set grid
set key top right
plot '2T-basic-Lee-intensity-at-peak.dat' using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  3 t '2T Basic Lee', \
     '2T-full-PLF-intensity-at-peak.dat'      using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  1 t '2T Full PLF', \
     '2T-full-PLFM-intensity-at-peak.dat'     using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  2 t '2T Full PLFM', \
     '3T-full-PLF-intensity-at-peak.dat'      using 1:($2*1.0e-4*1.0e-6) with lines lw 4 lt  4 t '3T Full PLF', \
     'E090806_IvW.txt'                    using 1:2 every ::230 w l lw 2 lt -1 t 'EAST Shot 46/18 UV'