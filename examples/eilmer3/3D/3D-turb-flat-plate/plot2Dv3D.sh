
gnuplot<<EOF
set term postscript eps enhanced 20
set output "mabey-cf1-2v3.eps"
set title "c_f profile along wall - 3D vs 2D"
set xlabel "x, m"
set ylabel "c_f"
set yrange [0:0.002]
set key left top
plot "./viscous_data0.dat" using (\$1):(\$5) title "Eilmer3-3D" with lines lt 1 lw 2, \
     "./viscous_data2D.dat" using (\$1):(\$5) title "Eilmer3-2D" with lines lt 2 lw 2, \
     "./experimental_data_cf.txt" using (\$1):(\$2):(0.10*0.001) \
     title "Mabey" with yerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20
set output "mabey-u-2v3.eps"
set title "Velocity profile near wall at x = 0.368 m - 3D vs 2D"
set xlabel "u-velocity, m/s"
set ylabel "z, mm"
set key left top
set xrange [0:800]
set yrange [0:25]
plot "./turb_flat_plate-x-368mm0.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-3D" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm2D.dat" using (\$6):((0.16-\$3)*1000) title "Eilmer3-2D" with lines lt 2 lw 2, \
     "./experimental_data.txt" using (\$7*712.89):((\$2)*1000):(0.05*712.89) \
     title "Mabey" with xerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20    
set output "mabey-temperature-2v3.eps"
set title "Temperature profile near wall at x = 0.368 m - 3D vs 2D"
set xlabel "Temperature, K"
set ylabel "z, mm"
set key right top
set xrange [0:350]
set yrange [0:25]
plot "./turb_flat_plate-x-368mm0.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-3D" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm2D.dat" using (\$20):((0.16-\$3)*1000) title "Eilmer3-2D" with lines lt 2 lw 2, \
     "./experimental_data.txt" using (\$8*62.157):((\$2)*1000):(0.05*62.157) \
     title "Mabey" with xerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20     
set output "mabey-pitot-2v3.eps"
set title "Pitot-pressure profile near wall at x = 0.368 m - 3D vs 2D"
set xlabel "P_{pitot}, kPa"
set ylabel "z, mm"
set key left top
set xrange [0:100]
set yrange [0:25]
plot "./turb_flat_plate-x-368mm0.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-3D" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm2D.dat" using ((\$22)/1000):((0.16-\$3)*1000) title "Eilmer3-2D" with lines lt 2 lw 2, \
     "./experimental_data.txt" using ((\$3)*3.1634):(\$2*1000):(0.05*84.57) \
     title "Mabey" with xerr pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20    
set output "mabey-mach-2v3.eps"
set title "Mach number profile near wall at x = 0.368 m - 3D vs 2D"
set xlabel "Mach number"
set ylabel "z, mm"
set key left top
set xrange [0:5]
set yrange [0:25]
plot "./turb_flat_plate-x-368mm0.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-3D" with lines lt 1 lw 2, \
     "./turb_flat_plate-x-368mm2D.dat" using (\$21):((0.16-\$3)*1000) title "Eilmer3-2D" with lines lt 2 lw 2, \
     "./experimental_data.txt" using ((\$6)*4.5099):((\$2)*1000):(0.05*4.5099) \
     title "Mabey" with xerr pt 6 lt 1
EOF
