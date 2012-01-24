# bw_post.sh

e3post.py --job=bw --tindx=9999 --vtk-xml --add-mach

# Plot the surface pressure on the wedge
# We want the EAST edge of block 1 and the SOUTH edge of block 1
e3post.py --job=bw --tindx=9999 --output-file=bw_surface.data \
    --slice-list="0,-1,:,0;1,:,0,0"
awk -f surface_pressure.awk bw_surface.data > bw_surface_p_coeff.data

gnuplot <<EOF
set term postscript eps enhanced 20
set output "bw_surface_pressure.eps"
set title "Blunted wedge: surface pressure coefficient."
set xlabel "s/R_n"
set ylabel "Pressure Coefficient, (p - p_{/Symbol \245})/q_{/Symbol \245}"
set yrange [0.0:2.0]
plot "bw_surface_p_coeff.data" using 1:2 title "CFD at t=399us" with lines, \
     "bw_surface_p_coeff.data" using 1:3 title "Modified Newtonian" with lines
EOF

# Plot the axial force coefficient.
awk -f xforce.awk bw.e3shared.log > bw_xforce.data

gnuplot <<EOF
set term postscript eps 20
set output "bw_xforce.eps"
set title "Blunted wedge: x-force history"
set xlabel "t, microseconds"
set ylabel "x-force, N"
set yrange [0:35000]
set key top left
plot "bw_xforce.data" using 1:2 title "total" with lines, \
     "bw_xforce.data" using 1:3 title "cylinder" with lines, \
     "bw_xforce.data" using 1:4 title "wedge" with lines
EOF


