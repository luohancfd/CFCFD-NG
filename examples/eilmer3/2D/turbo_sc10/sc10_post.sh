#! /bin/sh
# sc10_post.sh
# 2D sc10 profile, extract data and plot it.

# Around the blade
e3post.py --job=sc10 --output-file=surface.dat --tindx=9999 \
    --slice-list="0,:,0,0;1,-1,:,0;2,:,-1,0;3,:,-1,0"
#                 0south  1east    2north   3north

# Extract the solution data over whole flow domain and reformat.
e3post.py --job=sc10 --vtk-xml --add-mach

# Calculate average flow properties at inlet and outlet
turbo_post.py sc10

gnuplot <<EOF
set term postscript eps 20
set output "surface_p.eps"
set title "Standard Condition 10 M=0.7"
set xlabel "x, m"
set ylabel "Pressure ratio p/p0"
#set xrange [0.05:0.8]
#set yrange [-0.6:1.0]
plot "surface.dat" using 1:(\$9/100000) title "Elmer3" with points, \
     "rpmTurboSc10SubsonicSteady-M0.7.txt" using 2:4 title "RPM-Turbo" with points
EOF

echo At this point, we should have a plotted data.

