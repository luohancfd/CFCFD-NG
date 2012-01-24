# extract-history.awk
BEGIN {
    print "# t(ms)  Mach  p(kPa)";
}

$2==59 && $3==1 {
    t = $1;
    u = $10;
    v = $11;
    a = $14;
    p = $13;
    print t*1000.0, sqrt(u*u+v*v)/a, p/1000.0;
}
