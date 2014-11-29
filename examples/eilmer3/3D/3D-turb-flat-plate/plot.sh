gnuplot<<EOF
set term postscript eps enhanced 20
set output "mabey-cf-3D.eps"
set title "c_f profile along wall - 3D"
set xlabel "x, m"
set ylabel "c_f"
set yrange [0:0.002]
set key right top
plot "./viscous_data0.dat" using (\$1):(\$5) title "Eilmer3-cr0" with lines lt 1 lw 2, \
     "./viscous_data1.dat" using (\$1):(\$5) title "Eilmer3-cr1" with lines lt 1 lw 2, \
     "./viscous_data2.dat" using (\$1):(\$5) title "Eilmer3-cr2" with lines lt 1 lw 2, \
     "./viscous_data3.dat" using (\$1):(\$5) title "Eilmer3-cr3" with lines lt 1 lw 2, \
     "./viscous_data4.dat" using (\$1):(\$5) title "Eilmer3-cr4" with lines lt 1 lw 2, \
     "./viscous_data5.dat" using (\$1):(\$5) title "Eilmer3-cr5" with lines lt 1 lw 2, \
     "./viscous_data6.dat" using (\$1):(\$5) title "Eilmer3-cr6" with lines lt 1 lw 2, \
     "./viscous_data7.dat" using (\$1):(\$5) title "Eilmer3-cr7" with lines lt 1 lw 2, \
     "./viscous_data8.dat" using (\$1):(\$5) title "Eilmer3-cr8" with lines lt 1 lw 2, \
     "./viscous_data9.dat" using (\$1):(\$5) title "Eilmer3-cr9" with lines lt 1 lw 2, \
     "./experimental_data_cf.txt" using (\$1):(\$2):(0.10*0.001) \
     title "Mabey" with yerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20
set output "mabey-u-3D.eps"
set title "Velocity profile near wall at x = 0.368 m - 3D"
set xlabel "u-velocity, m/s"
set ylabel "z, mm"
set key left top
set xrange [0:800]
set yrange [0:25]
plot "./turb_flat_plate-x-368mm0.dat" using (\$6):((0.16-\$3)*1000) \
title "Eilmer3-cr0" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm1.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-cr1" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm2.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-cr2" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm3.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-cr3" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm4.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-cr4" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm5.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-cr5" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm6.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-cr6" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm7.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-cr7" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm8.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-cr8" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm9.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-cr9" with lines lt 1 lw 2, \
     "./experimental_data.txt" using (\$7*712.89):((\$2)*1000):(0.05*712.89) \
     title "Mabey" with xerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20    
set output "mabey-temperature-3D.eps"
set title "Temperature profile near wall at x = 0.368 m - 3D"
set xlabel "Temperature, K"
set ylabel "z, mm"
set key right top
set xrange [0:350]
set yrange [0:25]
plot "./turb_flat_plate-x-368mm0.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-cr0" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm1.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-cr1" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm2.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-cr2" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm3.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-cr3" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm4.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-cr4" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm5.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-cr5" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm6.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-cr6" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm7.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-cr7" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm8.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-cr8" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm9.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-cr9" with lines lt 1 lw 2, \
     "./experimental_data.txt" using (\$8*62.157):((\$2)*1000):(0.05*62.157) \
     title "Mabey" with xerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20     
set output "mabey-pitot-3D.eps"
set title "Pitot-pressure profile near wall at x = 0.368 m - 3D"
set xlabel "P_{pitot}, kPa"
set ylabel "z, mm"
set key left top
set xrange [0:100]
set yrange [0:25]
plot "./turb_flat_plate-x-368mm0.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-cr0" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm1.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-cr1" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm2.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-cr2" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm3.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-cr3" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm4.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-cr4" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm5.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-cr5" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm6.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-cr6" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm7.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-cr7" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm8.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-cr8" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm9.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-cr9" with lines lt 1 lw 2, \
     "./experimental_data.txt" using ((\$3)*3.1634):(\$2*1000):(0.05*84.57) \
     title "Mabey" with xerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20    
set output "mabey-mach-3D.eps"
set title "Mach number profile near wall at x = 0.368 m - 3D"
set xlabel "Mach number"
set ylabel "z, mm"
set key left top
set xrange [0:5]
set yrange [0:25]
plot "./turb_flat_plate-x-368mm0.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-cr0" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm1.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-cr1" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm2.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-cr2" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm3.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-cr3" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm4.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-cr4" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm5.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-cr5" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm6.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-cr6" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm7.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-cr7" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm8.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-cr8" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm9.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-cr9" with lines lt 1 lw 2, \
     "./experimental_data.txt" using ((\$6)*4.5099):((\$2)*1000):(0.05*4.5099) \
     title "Mabey" with xerr pt 6 lt 1
EOF
