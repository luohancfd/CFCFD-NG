BEGIN { print "# 1: t     2: L1 norm for rho" }
/Compare with reference function/ { t = $7 }
/rho L1/ { print t, $3 }

