#! /bin/sh
# sc10_new_post.sh
# 2D sc10_new profile, extract data and plot it.

# Around the blade
e3post.py --job=sc10_new --output-file=surface_new.dat --tindx=9999 \
    --slice-list="0,:,0,0;1,-1,:,0;2,:,-1,0;3,:,-1,0"
#                 0south  1east    2north   3north

# Extract the solution data over whole flow domain and reformat.
e3post.py --job=sc10_new --vtk-xml --add-mach --tindx=all

# Calculate average flow properties at inlet and outlet
turbo_post.py sc10_new

gnuplot <<EOF
set term postscript eps 20
set output "surface_p_new.eps"
set title "Standard Condition 10 M=0.7"
set xlabel "x, m"
set ylabel "Pressure ratio p/p0"
#set xrange [0.05:0.8]
#set yrange [-0.6:1.0]
plot "surface_new.dat" using 1:(\$9/100000) title "Eilmer3" with points, \
     "rpmTurboSc10SubsonicSteady-M0.7.txt" using 2:4 title "RPM-Turbo" with points
EOF

echo At this point, we should have a plotted data.

