#! /bin/sh
# cone20_run.sh
# exercise the Navier-Stokes solver for the cone20 test case.
# It is assumed that the path is set correctly.

# Prepare the simulation input files (parameter, grid and initial flow data).
# The SVG file provides us with a graphical check on the geometry
e3prep.py --job=cone20 --do-svg
if [ "$?" -ne "0" ] ; then
    echo "e3prep.py ended abnormally."
    exit
fi

# Integrate the solution in time, 
# recording the axial force on the cone surface.
time e3shared.exe -f cone20 --run --verbose
if [ "$?" -ne "0" ] ; then
    echo "e3shared.exe ended abnormally."
    exit
fi

# Extract the solution data and reformat.
# If no time is specified, the final solution found is output.
e3post.py --job=cone20 --vtk-xml

# Extract the average coefficient of pressure from the axial force
# records that were written to the simulation log file.
awk -f cp.awk e3shared.log > cone20_cp.dat

# Plot the average coefficient of pressure on the cone surface.
# We assume that the high-resolution data file is also available.
gnuplot <<EOF
set term postscript eps enhanced 20
set output "cone20_cp.ps"
set style line 1 linetype 1 linewidth 3.0 
set title "20 degree cone in Mach 1.5 flow"
set xlabel "time, ms"
set ylabel "average C_p"
set xtic 1.0
set ytic 0.1
set yrange [0:0.5]
set key bottom right
set arrow from 5.2,0.387 to 5.8,0.387 nohead linestyle 1
set label "Value from\nNACA 1135\nChart 6" at 5.0,0.3 right
set arrow from 5.0,0.3 to 5.5,0.387 head
plot "cone20_cp.dat" using 1:2 title "10x40+30x40", \
     "cone20_cp_hi-res.dat" using 1:2 title "20x80+60x80" with lines
EOF

echo "At this point, we should have a solution that can be viewed with Paraview."
