#! /bin/bash
# plot.sh
e3post.py --job=swlbli --tindx=last --add-mach --output-file=bl.data \
    --slice-list="2,:,0,0;4,:,0,0;6,:,0,0;8,:,0,0;10,:,0,0;12,:,0,0;14,:,0,0;16,:,0,0;18,:,0,0"

gnuplot <<EOF
set term postscript eps 20
set output "pressure.eps"
set title "Static pressure along the plate pf/p0=1.4"
set ylabel "p, kPa"
set yrange [0:10]
set xlabel "x, mm"
set key left top
plot "./bl.data" using (\$1*1000.0):(\$9/1000.0) title "Eilmer3" with lines, \
     "./notes/fig6b-pressure.data" using (\$1*25.4):(\$2*49.64) title "Hakkinen Fig.6b" with points pt 4
EOF
echo "At this point, we should have data to view"


