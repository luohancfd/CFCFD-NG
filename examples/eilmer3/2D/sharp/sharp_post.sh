#! /bin/sh
# sharp_post.sh
# Sharp 2D body, extract data and plot it.

# Extract the solution data over whole flow domain and reformat.
e3post.py --job=sharp --tindx=15 --vtk-xml

# Extract surface pressure and plot.
e3post.py --job=sharp --output-file=sharp_surface.dat --tindx=15 \
    --slice-list="1,:,0,0"

gnuplot <<EOF
set term postscript eps 20
set output "sharp_surface_p.eps"
set title "Sharp 2D Body in Mach 3 Freestream"
set xlabel "x, m"
set ylabel "Pressure, kPa"
set xrange [0.0:10.0]
set yrange [0.0:800]
plot "sharp_surface.dat" using 1:(\$9/1000) with lines
EOF

echo At this point, we should have a plotted data.

