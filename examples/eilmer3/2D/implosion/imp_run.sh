#!/bin/sh
# imp_run.sh
e3prep.py --job=imp
e3shared.exe --job=imp --run
e3post.py --job=imp --vtk-xml --tindx=all --add-mach

e3post.py --job=imp --tindx=9999 --add-mach --output-file=xaxis_profile.data \
    --slice-along-line="0.0,0.0,0.0,1.0,0.0,0.0,99"
e3post.py --job=imp --tindx=9999 --add-mach --output-file=diagonal_profile.data \
    --slice-along-line="0.0,0.0,0.0,1.0,1.0,0.0,100"
e3post.py --job=imp --tindx=9999 --add-mach --output-file=degree30_profile.data \
    --slice-along-line="0.0,0.0,0.0,0.8660,0.50,0.0,100"

gnuplot <<EOF
set term postscript eps 20
set output "density-vs-radius.eps"
set title "MNM Implosion Problem"
set xlabel "r/L"
set ylabel "Normalised density"
set xrange [0.0:1.0]
set yrange [0.0:12.0]
plot "xaxis_profile.data" using (sqrt(\$1*\$1+\$2*\$2+\$3*\$3)/1.0):(\$5/1.161) \
     title "x-axis" with lines, \
     "diagonal_profile.data" using (sqrt(\$1*\$1+\$2*\$2+\$3*\$3)/1.0):(\$5/1.161) \
     title "45 degrees" with lines, \
     "degree30_profile.data" using (sqrt(\$1*\$1+\$2*\$2+\$3*\$3)/1.0):(\$5/1.161) \
     title "30 degrees" with lines
EOF