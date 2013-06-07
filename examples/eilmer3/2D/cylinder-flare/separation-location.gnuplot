# separation-location.gnuplot

xzero(t) = xf + dx * exp(-t/tau)
xf = 59.85
dx = 46.56
tau = 0.563
# Couldn't get gnuplot to behave for fitting the function.
# Have used scipy instead.
# fit xzero(t) './separation-location.data' using ($1*1000):($2*1000) via xf, dx, tau

set term postscript eps 20
set output 'separation-location.eps'
set title 'Cylinder with extended flare, separation location'
set ylabel 'xzero, mm'
set xlabel 't, ms'
set key top right
plot './separation-location.data' using ($1*1000):($2*1000) \
    title 'Eilmer3' with points pt 4, \
    xzero(x) title '59.85+46.56*exp(-t/0.563)' with lines ls 1 lw 2 
