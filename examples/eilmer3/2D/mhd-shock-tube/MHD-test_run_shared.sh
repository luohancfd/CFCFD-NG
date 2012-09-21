#! /bin/sh
# MHD-test_run_mpi.sh
# exercise the Navier-Stokes solver with added MHD capability for the MHD-test test case.
# It is assumed that the path is set correctly.

# Prepare the simulation input files (parameter, grid and initial flow data).
# The SVG file provides us with a graphical check on the geometry
e3prep.py --job=MHD-test --do-svg
if [ "$?" -ne "0" ] ; then
    echo "e3prep.py ended abnormally."
    exit
fi

# Integrate the solution in time 
time e3shared.exe -f MHD-test --run --verbose
if [ "$?" -ne "0" ] ; then
    echo "e3shared.exe ended abnormally."
    exit
fi

# Extract the solution data and reformat.
# If no time is specified, the final solution found is output.
e3post.py --job=MHD-test --vtk-xml --tindx=all

echo "At this point, we should have a solution that can be viewed with Paraview."


# plot some stuff
# extract a line of data
#e3post.py --job=MHD-test --slice-along-line="-1.0,0.0,0.0,1.0,0.0,0.0,1000" --tindx=9999

gnuplot <<EOF
set term postscript eps enhanced 20
set output "MHD_shock_density.ps"
set style line 1 linetype 1 linewidth 3.0 
set title "MHD Shock Tube"
set xlabel "x, m"
set ylabel "density, kg/m^3"
set xtic 0.5
set ytic 0.1
set yrange [0.1:0.6]
set key top right
plot "exact-sol.dat" using 1:2 title 'rho-exact' with lines,\
     "profile.data" using 1:5 title 'rho-MHD' with lines
EOF

gnuplot <<EOF
set term postscript eps enhanced 20
set output "MHD_shock_B_y.ps"
set style line 1 linetype 1 linewidth 3.0 
set title "MHD Shock Tube"
set xlabel "x, m"
set ylabel "B.y"
set xtic 0.5
set ytic 50
set yrange [10.0:330.0]
set key bottom left
plot "exact-sol.dat" using 1:3 title 'B.y-exact' with lines,\
     "profile.data" using 1:10 title 'B.y-MHD' with lines
EOF

gnuplot <<EOF
set term postscript eps enhanced 20
set output "MHD_shock_B_z.ps"
set style line 1 linetype 1 linewidth 3.0 
set title "MHD Shock Tube"
set xlabel "x, m"
set ylabel "B.z"
set xtic 0.5
set ytic 50
set yrange [-10.0:375.0]
set key bottom right
plot "exact-sol.dat" using 1:4 title 'B.z-exact' with lines,\
     "profile.data" using 1:11 title 'B.z-MHD' with lines
EOF


