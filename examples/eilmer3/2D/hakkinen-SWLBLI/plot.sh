#! /bin/bash
# plot.sh
e3post.py --job=swlbli --tindx=last --add-mach --output-file=bl.data \
    --slice-list="2,:,0,0;4,:,0,0;6,:,0,0;8,:,0,0;10,:,0,0;12,:,0,0;14,:,0,0;16,:,0,0;18,:,0,0"

gnuplot pressure.gnuplot
gnuplot shear.gnuplot
echo "At this point, we should have pictures to view"


