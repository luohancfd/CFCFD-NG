# all basic plots

set term postscript eps enhanced color "Helvetica" 32
set size 1.75, 1.5
set output "profile_T-C.eps"
set xlabel "Distance behind shock, x (mm)"
set ylabel "Temperature (K)"
set xrange [0:20]
set yrange [0:13000]
set grid
set key bottom right
plot "TC2M1-2T-basic-Lee.data"      u ($1*1000.0):30 t '2T Basic Lee' with lines lt  2 lw 8, \
     "TC2M1-2T-full-Lee.data"       u ($1*1000.0):30 t '2T Full Lee' with lines lt  5 lw 8, \
     "TC2M1-2T-full-Park.data"      u ($1*1000.0):30 t '2T Full Park' with lines lt  8 lw 8, \
     "TC2M1-2T-full-PL.data"        u ($1*1000.0):30 t '2T Full PL'  with lines lt  4 lw 8, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):30 t '2T Full PLF' with lines lt  -1 lw 4, \
     "TC2M1-2T-full-PLFM.data"      u ($1*1000.0):30 t '2T Full PLFM' with lines lt  0 lw 8, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):30 t '3T Full PLF'   with lines lt  1 lw 8, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):43 t '3T Full PLF: T_{e}'   with lines lt 3 lw 8, \
     "E091206-spec_opt-results.txt" u ($1+1):3:5:(1.5*$6*$3) w xyerrorbars ps 2 pt 3 lt -1 t 'EAST 46-19: T_{v,CN}'
     
set output "mole-profiles.eps"
set ylabel "Species number density, n (particles/cm^3)"
set xlabel "distance, x (mm)"
set logscale y
unset logscale x
set yrange [1.0e12:1.0e18]
set xrange [-1:70.0]
set key bottom left
set grid 
set label "CO"   at 65, 1.4518e-1*6.022e17
set label "CN"   at 65, 1.5739e-4*6.022e17
set label "C_2"  at 65, 9.6879e-5*6.022e17
set label "e^-"  at 65, 3.9389e-3*6.022e17
plot "TC2M1-2T-full-Lee.data"  u ($1*1000.0):($10*6.022e17)  t "2T Full Lee"  w l lt  5 lw 4, \
     "TC2M1-2T-full-Park.data" u ($1*1000.0):($10*6.022e17)  t "2T Full Park" w l lt  8 lw 4, \
     "TC2M1-2T-full-PL.data"   u ($1*1000.0):($10*6.022e17)  t "2T Full PL"   w l lt  4 lw 4, \
     "TC2M1-2T-full-PLF.data"  u ($1*1000.0):($10*6.022e17)  t "2T Full PLF"  w l lt -1 lw 4, \
     "TC2M1-2T-full-PLFM.data" u ($1*1000.0):($10*6.022e17)  t "2T Full PLF"  w l lt  1 lw 4, \
     "CEA-mfs.txt"             u           1:($3*6.022e17)  notitle          w l lt  -1 lw 3, \
     "TC2M1-2T-full-Lee.data"  u ($1*1000.0):($12*6.022e17) notitle          w l lt  5 lw 4, \
     "TC2M1-2T-full-Park.data" u ($1*1000.0):($12*6.022e17) notitle          w l lt  8 lw 4, \
     "TC2M1-2T-full-PL.data"   u ($1*1000.0):($12*6.022e17) notitle          w l lt  4 lw 4, \
     "TC2M1-2T-full-PLF.data"  u ($1*1000.0):($12*6.022e17) notitle          w l lt -1 lw 4, \
     "TC2M1-2T-full-PLFM.data" u ($1*1000.0):($12*6.022e17) notitle          w l lt  1 lw 4, \
     "CEA-mfs.txt"             u           1:($5*6.022e17)  notitle          w l lt  -1 lw 3, \
     "TC2M1-2T-full-Lee.data"  u ($1*1000.0):($15*6.022e17) notitle          w l lt  5 lw 4, \
     "TC2M1-2T-full-Park.data" u ($1*1000.0):($15*6.022e17) notitle          w l lt  8 lw 4, \
     "TC2M1-2T-full-PL.data"   u ($1*1000.0):($15*6.022e17) notitle          w l lt  4 lw 4, \
     "TC2M1-2T-full-PLF.data"  u ($1*1000.0):($15*6.022e17) notitle          w l lt -1 lw 4, \
     "TC2M1-2T-full-PLFM.data" u ($1*1000.0):($15*6.022e17) notitle          w l lt  1 lw 4, \
     "CEA-mfs.txt"             u           1:($8*6.022e17)   notitle          w l lt  -1 lw 3, \
     "TC2M1-2T-full-Lee.data"  u ($1*1000.0):($29*6.022e17) notitle          w l lt  5 lw 4, \
     "TC2M1-2T-full-Park.data" u ($1*1000.0):($29*6.022e17)  notitle         w l lt  8 lw 4, \
     "TC2M1-2T-full-PL.data"   u ($1*1000.0):($29*6.022e17)  notitle         w l lt  4 lw 4, \
     "TC2M1-2T-full-PLF.data"  u ($1*1000.0):($29*6.022e17)  notitle         w l lt -1 lw 4, \
     "TC2M1-2T-full-PLFM.data" u ($1*1000.0):($29*6.022e17)  notitle         w l lt  1 lw 4, \
     "CEA-mfs.txt"             u           1:($22*6.022e17)  notitle         w l lt -1 lw 3

set term postscript eps enhanced color "Helvetica" 24
set size 2.0, 1.5
set output "profile_T-A.eps"
set xlabel "Distance behind shock, x (mm)"
set ylabel "Temperature (K)"
set xrange [-1:61.0]
set yrange [0:25000]
set y2label "Electron density (cm^{-3})"
set y2tics
set grid
plot "TC2M1-2T-basic-Lee.data"      u ($1*1000.0):2  t '2T Basic Lee: T'     with lines lt -1 lw 2, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):2  t '2T Full PLF: T'     with lines lt -1 lw 4, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):2  t '3T Full PLF: T'       with lines lt -1 lw 6, \
     "TC2M1-2T-basic-Lee.data"      u ($1*1000.0):30 t '2T Basic Lee: T_{v}' with lines lt  1 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):30 t '2T Full PLF: T_{v}' with lines lt  1 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):30 t '3T Full PLF: T_{v}'   with lines lt  1 lw 5, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):43 t '3T Full PLF: T_{e}'   with lines lt  2 lw 5, \
     "TC2M1-2T-basic-Lee.data"      u ($1*1000.0):($29*6.022e17) axes x1y2 t '2T Basic Lee: n_{e}' with lines lt  3 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($29*6.022e17) axes x1y2 t '2T Full PLF: n_{e}' with lines lt  3 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($29*6.022e17) axes x1y2 t '3T Full PLF: n_{e}'   with lines lt  3 lw 5
     
set term postscript eps enhanced color "Helvetica" 26
set size 1.75, 1.5
set output "profile_T-B.eps"
set xlabel "Distance behind shock, x (mm)"
set ylabel "Temperature (K)"
set xrange [0:20.0]
set yrange [0:25000]
set y2range [0:5.0e15]
set y2label "Electron number density, n_e (cm^{-3})"
set y2tics
set grid
# set key 11,29000
set key top right
plot "TC2M1-2T-basic-Lee.data"      u ($1*1000.0):2  t '2T Basic Lee: T'     with lines lt -1 lw 2, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):2  t '2T Full PLF: T'     with lines lt -1 lw 4, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):2  t '3T Full PLF: T'       with lines lt -1 lw 6, \
     "TC2M1-2T-basic-Lee.data"      u ($1*1000.0):30 t '2T Basic Lee: T_{v}' with lines lt  1 lw 2, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):30 t '2T Full PLF: T_{v}' with lines lt  1 lw 4, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):30 t '3T Full PLF: T_{v}'   with lines lt  1 lw 6, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):43 t '3T Full PLF: T_{e}'   with lines lt  2 lw 6, \
     "TC2M1-2T-basic-Lee.data"      u ($1*1000.0):($29*6.022e17) axes x1y2 t '2T Basic Lee: n_{e}' with lines lt  3 lw 2, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($29*6.022e17) axes x1y2 t '2T Full PLF: n_{e}' with lines lt  3 lw 4, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($29*6.022e17) axes x1y2 t '3T Full PLF: n_{e}'   with lines lt  3 lw 6, \
     "E091206-spec_opt-results.txt" u 1:4:5:(1.5*$6*$4) w xyerrorbars ps 2 pt 4 lt  1 t 'EAST 46-19: T_{r,CN}', \
     "E091206-spec_opt-results.txt" u 1:3:5:(1.5*$6*$3) w xyerrorbars ps 2 pt 3 lt -1 t 'EAST 46-19: T_{v,CN}'
     
unset y2tics
unset y2range
unset y2label
set output "J-profile-A.eps"
set ylabel "emission coefficient, J (W/m^{3}-sr)"
set logscale y
set yrange [1e4:1e9]
set xrange [-1:61.0]
set key top right
plot "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($44) t '2T Basic Lee: Total emission' w l lt -1 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($44) t '2T Full PLF: Total emission' w l lt -1 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($44) t '3T Full PLF:   Total emission' w l lt -1 lw 5, \
     "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($45) t '2T Basic Lee: CO 4+'          w l lt  1 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($45) t '2T Full PLF: CO 4+'          w l lt  1 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($45) t '3T Full PLF: CO 4+'            w l lt  1 lw 5, \
     "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($51) t '2T Basic Lee: CN Violet'      w l lt  2 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($51) t '2T Full PLF: CN Violet'      w l lt  2 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($51) t '3T Full PLF: CN Violet'        w l lt  2 lw 5, \
     "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($64) t '2T Basic Lee: C'        w l lt  3 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($64) t '2T Full PLF: C'        w l lt  3 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($64) t '3T Full PLF: C'          w l lt  3 lw 5, \
     "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($66) t '2T Basic Lee: O'        w l lt  4 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($66) t '2T Full PLF: O'        w l lt  4 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($66) t '3T Full PLF: O'          w l lt  4 lw 5
     
set output "J-profile-B.eps"
set ylabel "emission coefficient, J (W/m^{3}-sr)"
set logscale y
set yrange [1e4:1e9]
set xrange [-1:20.0]
set key top right
plot "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($44) t '2T Basic Lee: Total emission' w l lt -1 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($44) t '2T Full PLF: Total emission' w l lt -1 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($44) t '3T Full PLF:   Total emission' w l lt -1 lw 5, \
     "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($45) t '2T Basic Lee: CO 4+'          w l lt  1 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($45) t '2T Full PLF: CO 4+'          w l lt  1 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($45) t '3T Full PLF: CO 4+'            w l lt  1 lw 5, \
     "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($51) t '2T Basic Lee: CN Violet'      w l lt  2 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($51) t '2T Full PLF: CN Violet'      w l lt  2 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($51) t '3T Full PLF: CN Violet'        w l lt  2 lw 5, \
     "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($64) t '2T Basic Lee: C'        w l lt  3 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($64) t '2T Full PLF: C'        w l lt  3 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($64) t '3T Full PLF: C'          w l lt  3 lw 5, \
     "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($66) t '2T Basic Lee: O'        w l lt  4 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($66) t '2T Full PLF: O'        w l lt  4 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($66) t '3T Full PLF: O'          w l lt  4 lw 5
     
set output "J-profile-C.eps"
set ylabel "emission coefficient, J (W/m^{3}-sr)"
unset logscale y
set autoscale y
set xrange [-1:20.0]
set key top right
plot "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($51) t '2T Basic Lee: CN Violet' w l lt -1 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($51) t '2T Full PLF: CN Violet' w l lt -1 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($51) t '3T Full PLF: CN Violet' w l lt -1 lw 5, \
     "TC2M1-2T-basic-Lee.data"  u ($1*1000.0):($53) t '2T Basic Lee: C_2 Swan' w l lt 1 lw 1, \
     "TC2M1-2T-full-PLF.data"       u ($1*1000.0):($53) t '2T Full PLF: C_2 Swan' w l lt 1 lw 3, \
     "TC2M1-3T-full-PLF.data"       u ($1*1000.0):($53) t '3T Full PLF: C_2 Swan' w l lt 1 lw 5
     
set output "qrad-profiles.eps"
set ylabel "Half-range Radiative flux, q_{rad} (W/cm^2)"
unset logscale y
set yrange [0:2000]
set xrange [-1:61.0]
plot "2T-basic-Lee-qrad_V_x.dat" u ($1*1000):2 w lp lw 4 ps 3 t '2T Basic Lee', \
     "2T-full-PLF-qrad_V_x.dat" u ($1*1000):2 w lp lw 4 ps 3 t '2T Full PLF', \
     "2T-full-PLFM-qrad_V_x.dat" u ($1*1000):2 w lp lw 4 ps 3 t '2T Full PLFM', \
     "3T-full-PLF-qrad_V_x.dat" u ($1*1000):2 w lp lw 4 ps 3 t '3T Full PLFM'