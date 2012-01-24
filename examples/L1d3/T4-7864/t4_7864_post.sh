# t4_7864_post.sh
# We have now run the simulation further to include the
# primary diaphragm rupture and shock-processing of the
# test gas.

# (1) A space-time plot can be generated via a 
#     contour plotting package.
#     First, extract the pressure (on a logarithmic scale) and 
#     save it in a file (t4_7864_log_p.gen) suitable for 
#     the Generic contour  plotter.
#     Then contour variable 2 (0 is x, 1 is time) to get 
#     the X-T diagram.

sptime.exe -f t4_7864 -tstart 140.0e-3 -tstop 157.0e-3 -p -log -maxsol 1000
mb_cont.exe -fi t4_7864_log_p.gen -fo t4_7864_log_p.ps -ps \
            -var 2 -edge -notrueshape -xrange 24.0 37.0 2.0 \
            -yrange 143.0e-3 153.0e-3 2.0e-3

# (2) Flow history at specified locations can be extracted from
#     the history file t4_7864.hx and plotted with GNU-Plot.

l_hist.exe -f t4_7864 -tstop 157.0e-3 -xloc 0
l_hist.exe -f t4_7864 -tstop 157.0e-3 -xloc 1
l_hist.exe -f t4_7864 -tstop 157.0e-3 -xloc 2
l_hist.exe -f t4_7864 -tstop 157.0e-3 -xloc 3
l_hist.exe -f t4_7864 -tstop 157.0e-3 -xloc 4
l_hist.exe -f t4_7864 -tstop 157.0e-3 -xloc 5


